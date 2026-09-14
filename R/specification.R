#' GARCH Model Specification
#'
#' @description Specifies a GARCH model prior to estimation.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines.
#' @param y an xts vector with either \dQuote{Date} or \dQuote{POSIXct} index.
#' @param constant whether to estimate a constant (mean) for y,
#' @param model the type of GARCH model. Valid choices are \dQuote{garch} for
#' vanilla GARCH, \dQuote{gjrgarch} for asymmetric GARCH, \dQuote{egarch} for
#' exponential GARCH, \dQuote{aparch} for asymmetric power ARCH,
#' \dQuote{fgarch} for Family GARCH, \dQuote{cgarch} for component GARCH,
#' \dQuote{igarch} for integrated GARCH, and \dQuote{ewma} for the EWMA model.
#' @param order the (p,q) GARCH order.
#' @param arma the (ar,ma) order of the ARMA mean equation, jointly estimated
#' with the GARCH variance equation. Defaults to \sQuote{c(0,0)} (no ARMA
#' dynamics in the mean, matching prior behavior where the mean is either
#' zero or a constant). Supported for all models, including \dQuote{igarch}
#' and \dQuote{ewma}. Stationarity of the AR polynomial and invertibility
#' of the MA polynomial are guaranteed by construction via a Durbin-Levinson
#' (partial autocorrelation) reparameterization, and therefore do not
#' require additional nonlinear constraints during estimation. The estimated
#' AR/MA coefficients can be extracted with \code{\link{arma_coefficients}};
#' see also \code{\link{fitted}} for the (possibly time-varying) conditional
#' mean and \code{\link{residuals}} for the ARMA innovations. To fix an
#' entire AR and/or MA polynomial at target coefficients instead of
#' estimating them, set \sQuote{value} to the desired ar/ma coefficients
#' (not the internal pacf parameterization) and \sQuote{estimate = 0} on
#' every \sQuote{arpacf#}/\sQuote{mapacf#} row of the spec's
#' \sQuote{parmatrix} for that polynomial; the correct internal
#' Durbin-Levinson transform is then applied automatically. Fixing only
#' some (not all) of the lags of a given polynomial is not supported, since
#' the reparameterization couples all of a polynomial's lags together, and
#' will raise an error.
#' @param xreg an optional xts matrix of regressors in the conditional mean
#' equation, whose coefficients are named \sQuote{tau1..taum} in the parmatrix.
#' @param xreg_type the convention used for the regressor contribution:
#' \sQuote{arma_errors} (the default, matching \code{stats::arima}'s xreg
#' semantics) runs the AR/MA recursion on \sQuote{w_t = y_t - mu - x_t'tau} so
#' tau is the long-run marginal effect; \sQuote{armax} (the rugarch
#' convention) adds \sQuote{x_t'tau} to the conditional mean at time t only,
#' so tau is the impact effect. The two are algebraically identical whenever
#' the AR order is zero. Currently used by \code{estimate} only.
#' @param variance_targeting whether to use variance targeting rather than
#' estimating the conditional variance intercept.
#' @param vreg an optional xts matrix of regressors in the conditional variance
#' equation.
#' @param multiplicative whether to exponentiate the contribution of the
#' regressors else will be additive. In the case of the \dQuote{egarch} model,
#' since this is already a multiplicative model, the regressors are additive
#' irrespective of the choice made.
#' @param init the method to use to initialize the recursion of the conditional
#' variance.
#' @param backcast_lambda the decay power for the exponential smoothing used
#' when initializing the recursion using the backcast method.
#' @param sample_n the number of data points to use when initializing the
#' recursion using the sample method.
#' @param distribution a valid distribution from the available
#' re-parameterized distributions of the package.
#' @param ... not used.
#' @return An object of class \dQuote{tsgarch.spec}.
#' @aliases garch_modelspec
#' @rdname garch_modelspec
#' @author Alexios Galanos
#' @export
#'
#'
#'
garch_modelspec <- function(y, model = "garch", constant = FALSE,
                            order = c(1,1), arma = c(0,0), xreg = NULL,
                            xreg_type = c("arma_errors","armax"), variance_targeting = FALSE,
                            vreg = NULL, multiplicative = FALSE,
                            init = c("unconditional","sample","backcast"),
                            backcast_lambda = 0.7, sample_n = 10,
                            distribution = "norm", ...)

{
    # 1. check and initialize data
    parameter <- value <- NULL
    if  (!is.xts(y)) {
        stop("y must be an xts object")
    }
    if (NCOL(y) > 1) {
        stop("y must be a univariate time series")
    }
    series_name <- colnames(y)
    spec <- initialize_data(y)
    spec$target$series_name <- series_name
    # 2. validate arguments
    model <- match.arg(model[1], choices = valid_garch_models())
    distribution <- match.arg(distribution[1], choices = valid_distributions())
    init <- match.arg(init[1], choices = c("unconditional","sample","backcast"))
    multiplicative <- as.logical(multiplicative)
    constant <- as.logical(constant)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    if (sum(order) == 0) {
        variance_targeting <- FALSE
        init <- "unconditional"
    }
    # 2b. validate arma order
    if (length(arma) != 2 || any(!is.finite(arma)) || any(arma < 0) || any(arma != as.integer(arma))) {
        stop("arma must be a length 2 non-negative integer vector, e.g. c(1,1).")
    }
    arma <- as.integer(arma)

    # egarch already in logs
    if (model == "egarch") {
        if (multiplicative) warning("\nmultiplicative not valid for egarch model (already multiplicative due to log specification). Setting to FALSE")
        multiplicative <- FALSE
    }
    # cannot have variance targeting for igarch model
    if (model == "igarch") {
        if (variance_targeting) warning("\nvariance_targeting not possible in igarch model (Inf unconditional variance). Setting to FALSE")
        variance_targeting <- FALSE
    }
    if (model == "ewma") {
        if (variance_targeting) warning("\nvariance_targeting not valid for ewma (omega is zero). Setting to FALSE")
        variance_targeting <- FALSE
        multiplicative <- FALSE
    }
    order <- c(order[1], order[2])
    variance_targeting <- as.logical(variance_targeting)

    # 3. check regressors
    if (!is.null(vreg)) {
        vreg <- check_xreg(vreg, index(y))
        if (!multiplicative) {
            if (any(coredata(vreg) < 0)) {
                warning("\nvreg present with negative values and multiplcative = FALSE. Cannot guarantee positivity of variance.")
            }
        }
        if (variance_targeting) {
            warning("\nmultiplicative not available when variance_targeting = TRUE")
            multiplicative <- FALSE
        }
    } else {
        multiplicative <- FALSE
    }
    # 3b. check mean-equation regressors
    xreg_type <- match.arg(xreg_type[1], choices = c("arma_errors","armax"))
    if (!is.null(xreg)) {
        xreg <- check_xreg(xreg, index(y))
        # reject a rank-deficient design (including a constant column when the
        # constant is also estimated, which would be collinear with mu)
        xmat <- cbind(if (constant) rep(1, NROW(y)) else NULL, coredata(xreg))
        if (qr(xmat)$rank < NCOL(xmat)) stop("\nxreg is rank deficient (or, with constant = TRUE, contains a column collinear with the constant).")
    }
    # 4. populate specification object
    # cmodel: [maxpq, arch, garch, variance_targeting, multiplicative, distribution, ar, ma, armax]
    # maxpq is the overall pre-sample burn-in length, covering both the variance
    # recursion (garch order) and the mean recursion (arma order); the last flag
    # selects the armax (impact-effect) vs arma_errors (long-run) convention for
    # the mean-equation regressors.
    cmodel <- c(max(order, arma), order[1], order[2], as.integer(variance_targeting),
                as.integer(multiplicative), distribution_class(distribution),
                arma[1], arma[2], as.integer(xreg_type == "armax"))
    spec$model$model <- model
    # retain the user-facing model name separately since "ewma" is coerced to
    # "igarch" below, and re-specification helpers (e.g. tsbacktest) must be
    # able to recover the original choice
    spec$model$model_name <- model
    spec$model$order <- order
    spec$model$arma <- arma
    spec$model$variance_targeting <- variance_targeting
    spec$model$init <- init
    spec$model$backcast_lambda = backcast_lambda
    spec$model$sample_n <- sample_n
    if (is.null(vreg)) {
        spec$vreg$vreg <- matrix(0, ncol = 1, nrow = NROW(y))
        spec$vreg$include_vreg <- FALSE
        spec$vreg$multiplicative <- multiplicative
    } else {
        spec$vreg$vreg <- coredata(vreg)
        spec$vreg$include_vreg <- TRUE
        spec$vreg$multiplicative <- multiplicative
    }
    if (is.null(xreg)) {
        spec$xreg$xreg <- matrix(0, ncol = 1, nrow = NROW(y))
        spec$xreg$include_xreg <- FALSE
    } else {
        spec$xreg$xreg <- coredata(xreg)
        spec$xreg$include_xreg <- TRUE
    }
    spec$xreg$xreg_type <- xreg_type
    spec$distribution <- distribution
    # 5. populate parameters
    parmatrix <- initialize_parameters(model, y, constant = constant,
                                       order = order, arma = arma, xreg = xreg,
                                       variance_targeting = variance_targeting,
                                       vreg = vreg,
                                       multiplicative = multiplicative,
                                       init = init,
                                       backcast_lambda = backcast_lambda,
                                       sample_n = sample_n,
                                       distribution = distribution)
    if (model == "ewma") spec$model$model <- "igarch"
    if (sum(order) == 0) {
        parmatrix[parameter == "omega", value := as.numeric(var(y))]
    }
    spec$parmatrix <- parmatrix
    spec$model_options <- cmodel
    spec$model$constant <- constant
    class(spec) <- "tsgarch.spec"
    return(spec)
}


