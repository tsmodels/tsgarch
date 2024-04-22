#' Estimates an GARCH model given a specification object using maximum likelihood and autodiff
#'
#' @param object an object of class tsgarch.spec.
#' @param control solver control parameters.
#' @param solver only \dQuote{nloptr} is currently supported (see \code{\link[nloptr]{nloptr}}).
#' @param stationarity_constraint the bound on the inequality constraint for ensuring
#' the stationary of the GARCH process (see details).
#' @param ... not currently used.
#' @returns An object of class \dQuote{tsgarch.estimate}.
#' @details The underlying code is written using the TMB framework which uses
#' automatic differentiation and hence allows the generation of analytic
#' derivatives.
#' Stationarity is usually based on the condition that the persistence of the model
#' is less than 1. The argument \dQuote{stationarity_constraint} allows to fine tune
#' this. For example, setting it to a very high value will effectively render
#' this constraint inactive. The default of 0.999 has been found to be a reasonable
#' bound since values close to one may lead to problems.
#' Since the nloptr solver make use of analytic Jacobians for the inequality constraint,
#' these are either provided in closed form or calculated as part of the automatic
#' differentiation algorithms implemented in the package.
#' The estimation makes 2 passes to the solver. The first pass uses no parameter
#' scaling, whilst in the second pass the parameters (as well as bounds) are scaled
#' making use of the estimated hessian from the first pass in order to generate
#' a hopefully more robust solution.
#' @export estimate.tsgarch.spec
#' @aliases estimate
#' @method estimate tsgarch.spec
#' @rdname estimate
#' @author Alexios Galanos
#' @export
#'
estimate.tsgarch.spec <- function(object, solver = "nloptr", control = NULL, stationarity_constraint = 0.999, ...)
{
    # check for all fixed
    all_fixed_pars <- sum(object$parmatrix$estimate)
    if (all_fixed_pars == 0) {
        warning("\nall parameters are fixed (estimate = 0). Dispatching to tsfilter method instead.")
        out <- .filter.tsgarch.spec(object)
        return(out)
    }
    if (is.null(control)) control <- nloptr_fast_options(trace = FALSE)
    start_timer <- Sys.time()
    out <- .estimate_garch_model(object, solver, control, stationarity_constraint, ...)
    if (inherits(out, 'try-error')) {
        if (!object$model$variance_targeting) {
            warning("\ndefault model failed. Trying with variance targeting turned on.")
            new_spec <- .modelspec.tsgarch.spec(object, variance_targeting = TRUE)
            out <- estimate(new_spec, solver = solver, control = control, stationarity_constraint = stationarity_constraint, ...)
        }
    }
    end_timer <- Sys.time()
    elapsed <- end_timer - start_timer
    out$elapsed <- elapsed
    return(out)
}

#' Extract Model Coefficients
#'
#' @description Extract the estimated coefficients of a model.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @returns A numeric named vector of estimated coefficients.
#' @aliases coef
#' @method coef tsgarch.estimate
#' @rdname coef
#' @export
#'
#'
coef.tsgarch.estimate <- function(object, ...)
{
    out <- object$parmatrix[estimate == 1]$value
    names(out) <- object$parmatrix[estimate == 1]$parameter
    return(out)
}

#' Extract Volatility (Conditional Standard Deviation)
#'
#' @description Extract the conditional standard deviation from a GARCH model.
#' @param object an object of class \dQuote{tsgarch.estimate}, \dQuote{tsgarch.predict}
#' or \dQuote{tsgarch.simulate}.
#' @param ... not currently used.
#' @returns An xts vector of the conditional volatility.
#' @aliases sigma
#' @method sigma tsgarch.estimate
#' @rdname sigma
#' @export
#'
#'
sigma.tsgarch.estimate <- function(object, ...)
{
    out <- xts(object$sigma, object$spec$target$index)
    return(out)
}

#' Extract Model Fitted Values
#'
#' @description Extract the fitted values of the estimated model.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @returns An xts vector of the fitted values. Since only a constant is supported
#' in the conditional mean equation this is either a vector with a constant else
#' a vector with zeros.
#' @aliases fitted
#' @method fitted tsgarch.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.tsgarch.estimate <- function(object, ...)
{
    parameter <- NULL
    mu <- object$parmatrix[parameter == "mu"]$value
    idx <- object$spec$target$index
    f <- xts(rep(mu, length(idx)), idx)
    return(f)
}

#' Extract Model Residuals
#'
#' @description Extract the residuals of the estimated model.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param standardize logical. Whether to standardize the residuals by the
#' conditional volatility.
#' @param ... not currently used.
#' @returns An xts vector of the residuals. If the model had no constant in
#' the conditional mean equation then this just returns the original data (which
#' is assumed to be zero mean noise).
#' @aliases residuals
#' @method residuals tsgarch.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.tsgarch.estimate <- function(object, standardize = FALSE, ...)
{
    parameter <- NULL
    res <- object$spec$target$y_orig - object$parmatrix[parameter == "mu"]$value
    res <- xts(res, object$spec$target$index)
    if (standardize) {
        res <- res/sigma(object)
    }
    return(res)
}

#' The Covariance Matrix of the Estimated Parameters
#'
#' @param object an object of class tsgarch.estimate.
#' @param adjust logical. Should a finite sample adjustment be made? This amounts
#' to multiplication with n/(n-k) where n is the number of observations and k
#' the number of estimated parameters.
#' @param type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @returns The variance-covariance matrix of the estimated parameters.
#' @method vcov tsgarch.estimate
#' @aliases vcov
#' @rdname vcov
#' @export
#'
vcov.tsgarch.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    estimate <- NULL
    type <- match.arg(type[1],c("H","OP","QMLE","NW"))
    N <- nrow(estfun(object))
    if (type == "H") {
        V <- bread(object)
    } else if (type == "QMLE") {
        bread. <- bread(object)
        meat. <- meat_tsgarch(object, adjust = adjust)
        V <- bread. %*% meat. %*% bread.
    } else if (type == "OP") {
        V <- vcovOPG(object, adjust = adjust)
    } else if (type == "NW") {
        bread. <- bread(object)
        meat. <- meatHAC_tsgarch(object, adjust = adjust, ...)
        V <- bread. %*% meat. %*% bread.
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    colnames(V) <- rownames(V) <- par_names
    return(V)
}

#' Confidence Intervals for Model Parameters
#'
#' @param object an object of class tsgarch.estimate.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters
#' are considered.
#' @param level the confidence level required.
#' @param vcov_type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2
#' in % (by default 2.5% and 97.5%).
#' @method confint tsgarch.estimate
#' @aliases confint
#' @rdname confint
#' @export
#'
confint.tsgarch.estimate <- function(object, parm, level = 0.95, vcov_type = "H", ...)
{
    # extract the names of the estimated parameters
    estimate <- NULL
    par_names <- object$parmatrix[estimate == 1]$parameter
    coefficients <- coef(object)
    if (missing(parm)) {
        parm <- par_names
    } else if (is.numeric(parm)) {
        parm <- par_names[parm]
    } else if (is.character(parm)) {
        parm <- par_names[which(par_names %in% parm)]
    }
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    vc <- vcov(object, type = vcov_type)
    ses <- sqrt(diag(vcov(object, type = vcov_type)))
    ses <- ses[parm]
    ci[] <- coefficients[parm] + ses %o% fac
    return(ci)
}


#' Extract Log-Likelihood
#'
#' @description Extract the log likelihood of the model at the estimated optimal
#' parameter values.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return An object of class \dQuote{logLik} with attributes for \dQuote{nobs} and
#' \dQuote{df}. The latter is equal to the number of estimated parameters
#' plus 1 (the variance initialization value).
#' @aliases logLik
#' @method logLik tsgarch.estimate
#' @rdname logLik
#' @export
#'
#'
logLik.tsgarch.estimate <- function(object, ...)
{
    out <- -1.0 * object$loglik
    attr(out,"nobs") <- object$nobs
    attr(out,"df") <- object$npars
    class(out) <- "logLik"
    return(out)
}


#' GARCH Model Estimation Summary
#'
#' @description Summary method for class \dQuote{tsgarch.estimate}
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param vcov_type the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).
#' @param include_persistence whether to include the estimate of the persistence and
#' its calculated standard errors (calculated using the \code{\link[TMB]{sdreport}})
#' in the output.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @return A list with summary information of class \dQuote{summary.tsgarch.estimate}.
#' @aliases summary
#' @method summary tsgarch.estimate
#' @rdname summary
#' @export
#'
#'
summary.tsgarch.estimate <- function(object, digits = 4, vcov_type = "H", include_persistence = TRUE, ...)
{
    estimate <- NULL
    V <- vcov(object, type = vcov_type)
    est <- object$parmatrix[estimate == 1]$value
    par_names <- object$parmatrix[estimate == 1]$parameters
    se <- sqrt(diag(V))
    tval <- est/se
    coefficients <- cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = 2*(1 - pnorm(abs(tval))))
    persistence_summary <- t(data.frame("persistence" = object$persistence_summary))
    variance_target_summary <- t(data.frame("omega" = object$variance_target_summary, check.names = FALSE))
    if (object$spec$model$variance_targeting) {
        coefficients <- rbind(coefficients, variance_target_summary)
    }
    if (include_persistence) coefficients <- rbind(coefficients, persistence_summary)
    n_obs <- nobs(object)
    n_parameters <- length(coef(object))
    p <- persistence(object)
    llh <- -object$loglik
    elapsed <- object$elapsed
    conditions <- object$conditions[c("kkt1","kkt2","evratio")]
    distribution <- object$spec$distribution
    equation <- tsequation(object)
    coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
    uncvar <- unconditional(object)
    persist <- persistence(object)
    initvar <- object$var_initial
    setnames(coefficients, "rn","term")
    syms <- object$parmatrix[estimate == 1]$symbol
    if (object$spec$model$variance_targeting) {
        syms <- c(syms, "\\omega")
    }
    if (include_persistence) {
        syms <- c(syms,"P")
    }
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = llh, n_obs = n_obs, n_parameters = n_parameters,
                AIC = AIC(object),
                BIC = BIC(object),
                initvar = initvar,
                init_method = object$spec$model$init,
                variance_targeting = object$spec$model$variance_targeting,
                uncvar = uncvar,
                persistence = persist,
                include_persistence = include_persistence,
                elapsed = elapsed, conditions = conditions, equation = equation,
                model = object$spec$model$model, symbol = syms,
                equation = object$parmatrix[estimate == 1]$equation)
    class(out) <- "summary.tsgarch.estimate"
    return(out)
}

#' Model Estimation Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsgarch.estimate}
#' @param x an object of class \dQuote{summary.tsgarch.estimate}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param ... not currently used.
#' @return Invisibly returns the original summary object.
#' @aliases print.summary.tsgarch.estimate
#' @method print summary.tsgarch.estimate
#' @rdname print
#' @export
#'
#'
print.summary.tsgarch.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...)
{
    .print_screen(x, digits = digits, signif.stars = signif.stars, ...)
}

#' Transform a summary object into flextable
#' @description
#' Transforms a \dQuote{summary.tsgarch} object into a flextable
#' with options on symbolic representation and model equation.
#' @param x an object of class \dQuote{summary.tsgarch}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param include.symbols logical. If TRUE, replaces parameter names with their symbols (if they exist).
#' @param include.equation logical. If TRUE, adds a section with the symbolic model equation.
#' @param include.statistics logical. If TRUE, adds a section with summary statistics on the model.
#' @param table.caption an optional string for the table caption.
#' @param ... additional arguments passed to flextable method.
#' @importFrom flextable as_flextable
#' @return A flextable object.
#' @aliases as_flextable.summary.tsgarch.estimate
#' @method as_flextable summary.tsgarch.estimate
#' @rdname as_flextable.summary
#' @export
as_flextable.summary.tsgarch.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        signif.stars = getOption("show.signif.stars"),
                                        include.symbols = TRUE, include.equation = TRUE,
                                        include.statistics = TRUE,
                                        table.caption = paste0(toupper(x$model)," Model Summary"), ...)
{
    out <- .print_flextable(x, digits = digits, signif.stars = signif.stars,
                            include.symbols = include.symbols, include.equation = include.equation,
                            include.statistics = include.statistics,
                            table.caption = table.caption, ...)
    return(out)
}

#' Model Equation (LaTeX)
#'
#' @description Generates a list of model equations in LaTeX.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return A list of equations in LaTeX which can be used in documents. This is
#' a list with 4 slots for the conditional distribution, the conditional volatility,
#' the persistence and unconditional variance equations.
#' @details This method is called in the summary when the format output option
#' chosen is \dQuote{flextable}.
#' @aliases tsequation
#' @method tsequation tsgarch.estimate
#' @rdname tsequation
#' @export
#'
tsequation.tsgarch.estimate <- function(object, ...)
{
    if (object$spec$vreg$include_vreg) {
        vreg <- object$spec$vreg$vreg
    } else {
        vreg <- NULL
    }
    out <- switch(object$spec$model$model,
           "garch" =  .equation_garch(object$spec$model$order, vreg = vreg,
                                      multiplicative = object$spec$vreg$multiplicative,
                                      distribution = object$spec$distribution,
                                      variance_targeting = object$spec$model$variance_targeting),
           "egarch" =  .equation_egarch(object$spec$model$order, vreg = vreg,
                                      multiplicative = object$spec$vreg$multiplicative,
                                      distribution = object$spec$distribution,
                                      variance_targeting = object$spec$model$variance_targeting),
           "aparch" =  .equation_aparch(object$spec$model$order, vreg = vreg,
                                        multiplicative = object$spec$vreg$multiplicative,
                                        distribution = object$spec$distribution,
                                        variance_targeting = object$spec$model$variance_targeting),
           "gjrgarch" =  .equation_gjrgarch(object$spec$model$order, vreg = vreg,
                                        multiplicative = object$spec$vreg$multiplicative,
                                        distribution = object$spec$distribution,
                                        variance_targeting = object$spec$model$variance_targeting),
           "fgarch" =  .equation_fgarch(object$spec$model$order, vreg = vreg,
                                            multiplicative = object$spec$vreg$multiplicative,
                                            distribution = object$spec$distribution,
                                            variance_targeting = object$spec$model$variance_targeting),
           "cgarch" =  .equation_cgarch(object$spec$model$order, vreg = vreg,
                                        multiplicative = object$spec$vreg$multiplicative,
                                        distribution = object$spec$distribution,
                                        variance_targeting = object$spec$model$variance_targeting),
           "igarch" =  .equation_igarch(object$spec$model$order, vreg = vreg,
                                        multiplicative = object$spec$vreg$multiplicative,
                                        distribution = object$spec$distribution,
                                        variance_targeting = object$spec$model$variance_targeting))
    return(out)
}


#' Model Persistence
#'
#' @description General method the persistence of a model.
#' @param object an object of some class.
#' @param ... additional parameters passed to the method.
#' @return The persistence of the estimated model. For GARCH models, the formulation
#' varies by the type of model. See the vignette for more details.
#' @aliases persistence
#' @rdname persistence
#' @export
#'
#
persistence <- function(object, ...)
{
    UseMethod("persistence")
}

#' @method persistence tsgarch.estimate
#' @rdname persistence
#' @export
#'
#
persistence.tsgarch.estimate <- function(object, ...)
{
    env <- list()
    env$parmatrix <- copy(object$parmatrix)
    env$distribution <- object$spec$distribution
    pars <- object$parmatrix[estimate == 1]$value
    env$model <- object$spec$model$model
    return(.persistence(pars, env))
}

#' @method persistence tsgarch.spec
#' @rdname persistence
#' @export
#'
#
persistence.tsgarch.spec <- function(object, ...)
{
    env <- list()
    env$parmatrix <- copy(object$parmatrix)
    env$distribution <- object$distribution
    pars <- object$parmatrix[estimate == 1]$value
    env$model <- object$model$model
    return(.persistence(pars, env))
}

#' Unconditional Value
#'
#' @description Unconditional value of a GARCH model variance.
#' @param object an object of class \dQuote{tsgarch.estimate} or \dQuote{tsgarch.spec}.
#' @param ... not currently used.
#' @details
#' For some models, there is no closed form solution available for the unconditional
#' variance of higher order model (e.g. GARCH(2,1)) in which case a simulation
#' based approach is adopted to approximate the value.
#' @returns A numeric vector of length 1 of the unconditional variance of the model.
#' @aliases unconditional
#' @method unconditional tsgarch.estimate
#' @rdname unconditional
#' @export
#'
#
unconditional.tsgarch.estimate <- function(object, ...)
{
    out <- switch(object$spec$model$model,
           "garch" = .unconditional_garch(object),
           "egarch" = .unconditional_egarch(object),
           "aparch" = .unconditional_aparch(object),
           "gjrgarch" = .unconditional_gjrgarch(object),
           "fgarch" = .unconditional_fgarch(object),
           "cgarch" = .unconditional_cgarch(object),
           # technically Inf, but we calculate it for validation
           "igarch" = .unconditional_garch(object))
    return(out)
}


#' @method unconditional tsgarch.spec
#' @rdname unconditional
#' @export
#'
#
unconditional.tsgarch.spec <- function(object, ...)
{
    # trick to use same dispatch as the estimated object
    sp <- copy(object)
    sp$spec <- sp
    out <- switch(sp$spec$model$model,
                  "garch" = .unconditional_garch(sp),
                  "egarch" = .unconditional_egarch(sp),
                  "aparch" = .unconditional_aparch(sp),
                  "gjrgarch" = .unconditional_gjrgarch(sp),
                  "fgarch" = .unconditional_fgarch(sp),
                  "cgarch" = .unconditional_cgarch(sp),
                  # technically Inf, but we calculate it for validation
                  "igarch" = .unconditional_garch(sp))
    return(out)
}

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return A numeric value.
#' @aliases AIC
#' @method AIC tsgarch.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsgarch.estimate <- function(object, ..., k = 2)
{
    out <- ( -2.0 * as.numeric(logLik(object)) + k * object$npars)
    return(out)
}

#' Bayesian Information Criterion
#'
#' @description Extract the BIC from an estimated model.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases BIC
#' @method BIC tsgarch.estimate
#' @rdname BIC
#' @export
#'
#'
BIC.tsgarch.estimate <- function(object, ...)
{
    out <- -2 * as.numeric(logLik(object)) + object$npars * log(nobs(object))
    return(out)
}


#' Extract the Number of Observations
#'
#' @description Extract the number of observations from an estimated model.
#' This is principally intended to be used in computing BIC and used in other
#' tidy methods
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases nobs
#' @method nobs tsgarch.estimate
#' @rdname nobs
#' @export
#'
#'
nobs.tsgarch.estimate <- function(object, ...)
{
    return(object$nobs)
}


#' News Impact Curve
#'
#' @description General method the news impact of a model
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param epsilon a user supplied zero mean noise vector. If this is NULL
#' then a vector is created from the actual data using the minimum and maximum
#' range.
#' @param ... additional parameters passed to the method.
#' @note The method does not support higher order GARCH models.
#' @returns An object of class \dQuote{tsgarch.newsimpact}.
#' @aliases newsimpact
#' @rdname newsimpact
#' @export
#'
#
newsimpact <- function(object, epsilon = NULL, ...)
{
    UseMethod("newsimpact")
}

#' @method newsimpact tsgarch.estimate
#' @rdname newsimpact
#' @export
#'
#
newsimpact.tsgarch.estimate <- function(object, epsilon = NULL, ...)
{
    out <- switch(object$spec$model$model,
                  "garch" = .newsimpact_garch(object, epsilon),
                  "egarch" = .newsimpact_egarch(object, epsilon),
                  "aparch" = .newsimpact_aparch(object, epsilon),
                  "gjrgarch" = .newsimpact_gjrgarch(object, epsilon),
                  "fgarch" = .newsimpact_fgarch(object, epsilon),
                  "cgarch" = .newsimpact_cgarch(object, epsilon),
                  "igarch" = .newsimpact_igarch(object, epsilon))
    return(out)
}


#' News Impact Plot
#'
#' @description Plot method for newsimpact class.
#' @param x an object of class \dQuote{tsgarch.newsimpact}.
#' @param y not used.
#' @param ... additional arguments pass to \code{\link[graphics]{plot.xy}} other
#' than \dQuote{xlab}, \dQuote{ylab} and \dQuote{main}.
#' @method plot tsgarch.newsimpact
#' @rdname plot
#' @export
#'
#
plot.tsgarch.newsimpact <- function(x, y = NULL, ...)
{
    plot(x$x, x$y, ylab = x$yexpr, xlab = x$xexpr, main = paste0("News Impact [",x$model,"]"), ...)
    return(invisible(x))
}


#' Estimated Model Plots
#'
#' @description Plot method for \dQuote{tsgarch.estimate} class.
#' @param x an object of class \dQuote{tsgarch.estimate}.
#' @param y not used.
#' @param ... not used.
#' @method plot tsgarch.estimate
#' @rdname plot.tsgarch.estimate
#' @export
#'
#
plot.tsgarch.estimate <- function(x, y = NULL, ...)
{
    op <- par(no.readonly = TRUE)
    par(mgp = c(2,1,0))
    dpars <- extract_model_values(x, object_type = "estimate", "distribution")
    distribution <- x$spec$distribution
    dist_print <- distribution_abb(distribution)
    ylim <- c(pmin(min(abs(x$spec$target$y_orig)), min(as.numeric(sigma(x)))) + 1e-10,
              pmax(max(abs(x$spec$target$y_orig)), max(as.numeric(sigma(x)))))
    layout(matrix(c(1,1,2,3), ncol = 2, nrow = 2, byrow = T))
    par(mar = c(2.5, 4, 2, 3))
    plot(as.zoo(sigma(x)), ylab = expression(~sigma[t]), col = "steelblue", xlab = "", ylim = ylim,
         main = paste0(toupper(x$spec$model$model),"(",x$spec$model$order[1],",",x$spec$model$order[2],") - ",toupper(distribution)),
         cex.main = 0.9)
    lines(as.zoo(abs(x$spec$target$y)), col = "snow2")
    lines(as.zoo(sigma(x)), col = "steelblue")
    grid()
    par(mar = c(4, 4, 2, 3))
    plot(newsimpact(x), type = "l", col = 'gray5', cex.main = 0.9)
    grid()
    qqplot(as.numeric(residuals(x, standardize = T)), rdist(distribution = distribution, n = 30000, mu = 0, sigma = 1, skew = dpars[1], shape = dpars[2], lambda = dpars[3]),
           ylab = substitute(paste(nn, (z[t])), list(nn = dist_print)), xlab = expression(z[t] == epsilon[t]/sigma[t]),
           main = paste0("Standardized Residuals (z)\n Sampling Distribution : ",dist_print), col = "snow3", cex.main = 0.8)
    qqline(as.numeric(residuals(x, standardize = T)), col = "gray5", lty = 2)
    grid()
    par(op)
    return(invisible(x))
}


#' Model Filtering
#'
#' @description Filters new data based on an already estimated model or filters data
#' based on a specification object.
#' @param object an object of class \dQuote{tsgarch.estimate} or \dQuote{tsgarch.spec}.
#' @param y an xts vector of new values to filter. Can also be NULL in which case the
#` original object is returned (if of class \dQuote{tsgarch.estimate}), or the existing
#` data filtered (if of class \dQuote{tsgarch.spec}). See details.
#' @param newxreg not currently used,
#' @param newvreg variance regressors with the same number of rows as y. This can be either
#' a numeric or xts matrix. Only needed if the model was estimated with regressors in the
#' variance equation.
#' @param ... additional arguments for future expansion.
#' @return A \dQuote{tsgarch.estimate} object with updated information.
#' @details The method filters new data and updates the object with this new information
#' so that it can be called recursively as new data arrives. It is also possible to use
#' a specification object with fixed parameters, by appropriately setting the values
#' of the \dQuote{parmatrix} object in the specification slot. In this case, the returned object
#' will also be of classs \dQuote{tsgarch.estimate}.
#' If an object of \dQuote{tsgarch.spec} is used with y not NULL, then the method will
#' first filter the values of the data in the object, generating an object of
#' \dQuote{tsgarch.estimate} and then call the method again on this new object and the
#' new y values (and optionally any newvreg values). In this way, using either object
#' classes will return the exact same results. The timestamp indices of y must be
#' strictly greater than the maximum timestamp index of the data within the object (i.e.
#' we only filter on new data).
#' @aliases tsfilter
#' @method tsfilter tsgarch.estimate
#' @rdname tsfilter
#' @export
#'
#'
tsfilter.tsgarch.estimate <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, ...)
{
    return(.filter.tsgarch.estimate(object, y = y, newxreg = newxreg, newvreg = newvreg))
}

#' @aliases tsfilter
#' @method tsfilter tsgarch.spec
#' @rdname tsfilter
#' @export
#'
#'
tsfilter.tsgarch.spec <- function(object, y = NULL, newxreg = NULL, newvreg = NULL, ...)
{
    return(.filter.tsgarch.spec(object, y = y, newxreg = newxreg, newvreg = newvreg))
}


#' Model Prediction
#'
#' @description Prediction function for class \dQuote{tsgarch.estimate}.
#' @param object an object of class \dQuote{tsgarch.estimate}.
#' @param h the forecast horizon.
#' @param newxreg not currently used,
#' @param newvreg variance regressors rows equal to h. This can be either
#' a numeric or xts matrix. Only needed if the model was estimated with regressors in the
#' variance equation.
#' @param nsim the number of simulations to use for generating the simulated
#' predictive distribution. Defaults to zero (no simulated distribution).
#' @param sim_method the simulation method to use when \strong{nsim} great than zero. The \dQuote{parametric}
#' method samples from the model distribution whilst the \dQuote{bootstrap} from the standardized
#' model residuals.
#' @param block for the \dQuote{bootstrap} \strong{sim_method}, this allows to generate
#' block length samples (defaults to 1).
#' @param forc_dates an optional vector of forecast dates equal to h. If NULL will use the
#' implied periodicity of the data to generate a regular sequence of dates after the
#' last available date in the data.
#' @param init_states an optional vector of states to initialize the forecast.
#' If NULL, will use the last available state from the estimated model. This must
#' be equal to the max of the ARCH and GARCH terms.
#' @param seed an integer that will be used in a call to set.seed before simulating.
#' @param ... additional arguments for future expansion options.
#' @return A \dQuote{tsgarch.predict} object.
#' @details The bootstrap method considered here, is based on re-sampling innovations
#' from the empirical distribution of the fitted GARCH model to generate future
#' realizations of the series and sigma. This only considers distributional uncertainty
#' and will not generate prediction intervals for the 1-step ahead sigma forecast
#' for which only the parameter uncertainty is relevant in GARCH type models (and
#' not currently implemented).
#' When the horizon \strong{h} is equal to 1, no simulation is performaed since there is
#' no uncertainty to account for.
#' @references
#' \insertRef{Pascual2006}{tsgarch}
#' @aliases predict
#' @method predict tsgarch.estimate
#' @rdname predict
#' @export
#'
#'
predict.tsgarch.estimate <- function(object, h = 1, newxreg = NULL, newvreg = NULL, nsim = 0, sim_method = c("parametric","bootstrap"), block = 1, forc_dates = NULL, init_states = NULL, seed = NULL, ...)
{
    if (is.null(forc_dates)) {
        forc_dates <- .forecast_dates(forc_dates, h = h, sampling = object$spec$target$sampling, last_index = tail(object$spec$target$index, 1))
    }
    sim_method <- match.arg(sim_method[1], c("parametric","bootstrap"))
    h <- max(1, as.integer(h[1]))
    if (h == 1) nsim <- 0
    nsim <- max(0, as.integer(nsim[1]))
    block <- max(0, as.integer(block[1]))
    model <- object$spec$model$model
    p <- switch(model,
           "garch" = .predict_garch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "egarch" = .predict_egarch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "aparch" = .predict_aparch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "gjrgarch" = .predict_gjrgarch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "fgarch" = .predict_fgarch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "cgarch" = .predict_cgarch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...),
           "igarch" = .predict_garch(object = object, h = h, newxreg = newxreg, newvreg = newvreg, nsim = nsim, sim_method = sim_method, block = block, forc_dates = forc_dates, init_states = init_states, seed = seed, ...)
    )
    return(p)
}


#' Probability Integral Transform (PIT)
#'
#' @description Calculates and returns the conditional probability integral
#' transform given the data and estimated density
#' @param object an object.
#' @param ... not currently used.
#' @return An xts vector of the conditional probabilities.
#' @details The PIT is essentially the probabilities returned from the cumulative
#' distribution function (*p) given the data and estimated value of the mean,
#' conditional standard deviation and any other distributional parameters.
#' @aliases pit
#' @method pit tsgarch.estimate
#' @rdname pit
#' @export
#'
#
pit.tsgarch.estimate <- function(object, ...)
{
    parameter <- NULL
    distribution <- object$spec$distribution
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    sigma <- as.numeric(sigma(object))
    mu <- rep(object$parmatrix[parameter == "mu"]$value, length(sigma))
    r <- object$spec$target$y_orig
    p <- pdist(distribution, q = r, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda)
    p <- xts(p, object$spec$target$index)
    return(p)
}


#' Half Life
#'
#' @description Calculates and returns the half-life of a model.
#' @param object an object.
#' @param ... not currently used.
#' @details The half life is defined as the period it
#' takes a series  to reach half its long-term average values. For a GARCH model
#' this is defined as \eqn{log(0.5)/log(P)} where P is the persistence.
#' @return a numeric value representing the half life in periods based on the
#' frequency of the underlying data.
#' @aliases halflife
#' @method halflife tsgarch.estimate
#' @rdname halflife
#' @export
#'
#
halflife.tsgarch.estimate <- function(object, ...)
{
    p <- persistence(object)
    out <- abs(log(0.5)/log(p))
    return(out)
}

#' Omega (Variance Equation Intercept)
#'
#' @description Returns the intercept of a GARCH model.
#' @param object an object.
#' @param ... additional parameters passed to the method.
#' @details The intercept is either estimated directly as part of the model else
#' indirectly if variance targeting was selected.
#' @return a numeric value representing the value of the intercept.
#' @aliases omega
#' @rdname omega
#' @export
#'
#
omega <- function(object, ...)
{
    UseMethod("omega")
}

#' @method omega tsgarch.estimate
#' @rdname omega
#' @export
#'
#
omega.tsgarch.estimate <- function(object, ...)
{
    return(object$target_omega)
}

#' @method omega tsgarch.spec
#' @rdname omega
#' @export
#'
#
omega.tsgarch.spec <- function(object, ...)
{
    parameter <- NULL
    if (object$model$variance_targeting) {
        target_omega <- variance_target(object, ...)
    } else {
        target_omega <- object$parmatrix[parameter == "omega"]$value
    }
    return(target_omega)
}

#' Model Simulation
#'
#' @description Simulates paths of a GARCH model.
#' @param object an object of class \dQuote{tsgarch.spec}.
#' @param nsim the number of sample paths to generate.
#' @param seed an integer that will be used in a call to set.seed before simulating.
#' @param h the number of time steps to simulate paths for.
#' @param var_init the seed value for initializing the variance equation recursion.
#' If NULL, the variance target value is used based on the supplied parameters.
#' This should be a vector and assumes all sample paths are seeded the same way.
#' @param innov an optional matrix of dimensions nsim by h of zero mean unit variance
#' (standardized) innovations which will be used instead of the model distribution
#' for simulation. No checks are performed on whether the supplied values are
#' standardized.
#' @param innov_init an optional vector of initialization values for the
#' standardized innovations. This allows the seeding of the initial innovations
#' with user supplied values (useful when simulating forward from an existing
#' model for the purpose of continuing the modeled series from some fixed point).
#' @param vreg an optional vector of length h representing any pre-multiplied
#' variance regressors to use in the simulation.
#' @param burn burn in. Will be discarded before returning the output.
#' @param ... for aparch, fgarch, egarch and gjrgarch models, an optional
#' vector of length max(q,p) with values for initializing the ARCH equation and
#' named \dQuote{arch_initial}. This is mostly used for validation purposes. The
#' \dQuote{arch_initial} value is always returned by an estimated object.
#' @details Once a GARCH model is specified via \code{\link{garch_modelspec}},
#' the slot \dQuote{parmatrix} contains initial values for the parameters
#' which can be used to set them to any value for the simulation. This matrix
#' as well as details of the model (type, order, distribution) are the only
#' pieces of information used in the simulation. The \dQuote{vreg} argument in
#' the spec will be ignored. Instead, the user can supply a pre-multiplied vector
#' to the simulate function which will be used. Note that the \dQuote{multiplicative}
#' argument in the specification will be used in this case to determine how the
#' regressors enter the conditional variance equation. While the \dQuote{innov}
#' argument must be a matrix, all other values are vectors and assume that they
#' will be the same across all sample paths. If the user wants to assign different values
#' for arguments \dQuote{var_init}, \dQuote{innov_init} and \dQuote{vreg}, then
#' the simulate method should be called multiple times.
#' @return An object of class \dQuote{tsgarch.simulate} with slots for the
#' simulated sigma and series simulated distributions which are each of class
#' \dQuote{tsmodel.distribution}. The simulated error (not returned) is equal to
#' the simulated series less the mean equation constant if not equal to zero.
#' @aliases simulate
#' @method simulate tsgarch.spec
#' @rdname simulate
#' @export
#'
#'
simulate.tsgarch.spec <- function(object, nsim = 1, seed  = NULL, h = 1000, var_init = NULL,
                                  innov = NULL, innov_init = NULL, vreg = NULL, burn = 0, ...)
{
    new_spec <- .spec2newspec(object)
    new_spec$parmatrix$value <- object$parmatrix$value

    out <- switch(new_spec$model$model,
                  "garch" = .simulate_garch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                            innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "egarch" = .simulate_egarch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                              innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "aparch" = .simulate_aparch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                              innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "gjrgarch" = .simulate_gjrgarch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                              innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "fgarch" = .simulate_fgarch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                                  innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "cgarch" = .simulate_cgarch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                              innov_init = innov_init, vreg = vreg, burn = burn, ...),
                  "igarch" = .simulate_igarch(new_spec, h = h, seed = seed, nsim = nsim, var_init = var_init, innov = innov,
                                              innov_init = innov_init, vreg = vreg, burn = burn, ...))

    class(out) <- "tsgarch.simulate"
    return(out)
}
