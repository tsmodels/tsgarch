#' Model Parameter Profiling
#'
#' @description Profiles the model parameters under the assumptions of the model.
#' @details The function profiles the parameters of a model by simulating and
#' then estimating multiple paths from the assumed DGP. This makes it possible
#' to obtain a better understanding of the convergence properties (RMSE) of each
#' parameter under different data sizes.
#' @param object an object of class \dQuote{tsgarch.spec} with pre-set parameters.
#' @param nsim the number of paths to generate.
#' @param sizes a vector of data sizes for which to simulate and estimate.
#' @param var_init the variance to use to initialize the simulation.
#' @param seed an object specifying if and how the random number generator
#' should be initialized. See the simulate documentation for more details.
#' @param burn burn in samples.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... not currently used.
#' @note The function can use parallel functionality as long as the user has set
#' up a \code{\link[future]{plan}} using the future package.
#' External regressors are not supported at this time and an error will occur
#' if persent in the specification.
#' @returns An object of class \dQuote{tsgarch.profile}.
#' @aliases tsprofile
#' @method tsprofile tsgarch.spec
#' @rdname tsprofile
#' @export
#'
#'
tsprofile.tsgarch.spec <- function(object, nsim = 100, sizes = c(800, 1000, 1500, 2000, 3000), var_init = NULL, seed = NULL, burn = 0, trace = FALSE, ...)
{
    .error <- .squared_error <- .absolute_percent_error <- .absolute_error <- value <- NULL
    ":=" <- NULL
    if (object$vreg$include_vreg) {
        stop("\nexternal regressors in the variance equation (vreg) are not supported.")
    }
    value <- parameter <- actual <- NULL
    spec <- .spec2newspec(object)
    spec$parmatrix$value <- object$parmatrix$value
    object <- copy(spec)
    if (is.null(var_init)) {
        var_init = unconditional(object)
    }
    par_names <- object$parmatrix[estimate == 1]$parameter
    true_pars <- object$parmatrix[estimate == 1]$value
    true_persist <- persistence(object)
    true_skewness <- dskewness(object$distribution, skew = object$parmatrix[parameter == "skew"]$value,
                          shape = object$parmatrix[parameter == "shape"]$value,
                          lambda = object$parmatrix[parameter == "lambda"]$value)
    true_kurtosis <- dkurtosis(object$distribution, skew = object$parmatrix[parameter == "skew"]$value,
                          shape = object$parmatrix[parameter == "shape"]$value,
                          lambda = object$parmatrix[parameter == "lambda"]$value)
    # dkurtosis returns the excess kurtosis...add back 3
    true_kurtosis <- true_kurtosis + 3
    true_pars <- c(true_pars, true_persist, true_skewness, true_kurtosis)
    par_names <- c(par_names, "persistence","skewness","kurtosis")
    true_table <- data.table(parameter = par_names, actual = true_pars)
    sim <- simulate(object, nsim = nsim, h = max(sizes), var_init = var_init, seed = seed, burn = burn)
    nsim <- nrow(sim$series)
    sim_out <- list()
    if (trace) {
        prog_trace <- progressor(length(sizes))
    }
    start_time <- Sys.time()
    size_results <- lapply(1:length(sizes), function(j){
        if (trace) prog_trace()
        sim_out %<-% future_lapply(1:nsim, function(i){
            # set the threads to 1 for each fork
            data.table::setDTthreads(1)
            if (object$vreg$include_vreg) {
                v <- object$vreg$vreg[1:sizes[j],,drop = FALSE]
            } else {
                v <- NULL
            }
            s <- xts(sim$series[i,1:sizes[j]], as.Date(1:sizes[j]))
            spec_new <- garch_modelspec(s, model = object$model$model, distribution  = object$distribution,
                                        order = object$model$order, constant = object$model$constant,
                                        variance_targeting = object$model$variance_targeting,
                                        init = object$model$init, backcast_lambda = object$model$backcast_lambda,
                                        sample_n = object$model$sample_n, vreg = v, multiplicative = object$vreg$multiplicative)
            # account for any fixed parameters
            if (any(object$parmatrix$estimate == 0)) {
                fixed_pars <- object$parmatrix[estimate == 0]$parameter
                fixed_pars_value <- object$parmatrix[estimate == 0]$value
                spec_new$parmatrix[parameter %in% fixed_pars, value := fixed_pars_value]
                spec_new$parmatrix[parameter %in% fixed_pars, estimate := 0]
            }
            mod <- try(estimate(spec_new, control = nloptr_fast_options(trace = 0)), silent = TRUE)
            if (inherits(mod, 'try-error')) {
                coeff <- spec_new$parmatrix[estimate == 1]$value * as.numeric(NA)
                skewness <- kurtosis <- persist <- as.numeric(NA)
            } else{
                # check results and catch errors
                coeff <- unname(coef(mod))
                skewness <- dskewness(object$distribution, skew = mod$parmatrix[parameter == "skew"]$value,
                                      shape = mod$parmatrix[parameter == "shape"]$value,
                                      lambda = mod$parmatrix[parameter == "lambda"]$value)
                kurtosis <- dkurtosis(object$distribution, skew = mod$parmatrix[parameter == "skew"]$value,
                                      shape = mod$parmatrix[parameter == "shape"]$value,
                                      lambda = mod$parmatrix[parameter == "lambda"]$value)
                # dkurtosis returns the excess kurtosis...add back 3
                kurtosis <- kurtosis + 3
                persist <- mod$persistence_summary[1]
            }
            par_names <- spec_new$parmatrix[estimate == 1]$parameter
            return(data.table(parameter = c(par_names, "persistence", "skewness", "kurtosis"), value = c(coeff, persist, skewness, kurtosis),
                              draw = i, size = sizes[j]))
        },future.packages = c("tsmethods", "tsgarch", "xts", "data.table"))
        sim_out <- eval(sim_out)
        sim_out <- rbindlist(sim_out)
        return(sim_out)
    })
    size_results <- rbindlist(size_results)
    size_results <- merge(size_results, true_table, by = "parameter", all.x = TRUE)
    size_results[, .error := actual - value]
    size_results[, .squared_error := .error^2]
    size_results[, .absolute_error := abs(.error)]
    size_results[, .absolute_percent_error := .absolute_error/abs(actual)]
    summary_table <- size_results[,list(RMSE = sqrt(mean(.squared_error, na.rm = TRUE)), MAE = mean(.absolute_error, na.rm = TRUE),
                                        MAPE = mean(.absolute_percent_error, na.rm = TRUE),
                                        MEAN = mean(value, na.rm = TRUE),
                                        MEDIAN = mean(value, na.rm = TRUE),
                                        P20 = quantile(na.omit(value), 0.2),
                                        P80 = quantile(na.omit(value), 0.2)), by = c('parameter','size')]
    sol <- list()
    sol$profile <- size_results
    sol$summary <- summary_table
    sol$elapsed <- Sys.time() - start_time
    sol$model <- list(model = object$model$model, order = object$model$order, distribution = object$distribution,
                      parameters = par_names)
    class(sol) <- 'tsgarch.profile'
    return(sol)
}


#' GARCH Profile Summary
#'
#' @description Summary method for class \dQuote{tsgarch.profile}
#' @param object an object of class \dQuote{tsgarch.profile}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param measure either one of the included measure in the summary slot of
#' the returned object, which currently includes the relative error measures
#' \dQuote{RMSE}, \dQuote{MAE}, \dQuote{MAPE}, summary measures on the estimated
#' values \dQuote{MEAN}, \dQuote{MEDIAN}, \dQuote{P20} and \dQuote{P80},
#' else any other user calculated measure which has been generated in the summary table
#' post processing.
#' @param ... not currently used.
#' @return A list with summary information of class \dQuote{summary.tsgarch.profile}.
#' @aliases summary.tsgarch.profile
#' @method summary tsgarch.profile
#' @rdname summary.tsgarch.profile
#' @export
#'
#'
summary.tsgarch.profile <- function(object, digits = 4, measure = "RMSE", ...)
{
    tmp <- copy(object$summary)
    valid_measures <- colnames(tmp)[!colnames(tmp) %in% c("parameter","size")]
    if (!measure %in% valid_measures) stop("\nmeasure not in summary table!")
    actual <- object$profile[,list(actual = mean(actual)), by = c("parameter","size")]
    tmp <- merge(tmp, actual, by = c("parameter","size"))
    tab <- dcast(tmp, parameter + actual ~ size, value.var = measure)
    # re-order rows
    tab <- tab[object$model$parameters]
    colnames(tab)[1] <- "parameter/size"
    out <- list(table = tab, model = object$model, measure = measure)
    class(out) <- "summary.tsgarch.profile"
    return(out)
}

#' Profile Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsgarch.profile}
#' @param x an object of class \dQuote{summary.tsgarch.estimate.profile}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @return Invisibly returns the original summary object.
#' @aliases print.summary.tsgarch.profile
#' @method print summary.tsgarch.profile
#' @rdname print.summary.tsgarch.profile
#' @export
#'
#'
print.summary.tsgarch.profile <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    model <- paste0(toupper(x$model$model),"(",x$model$order[1],",",x$model$order[2],") - ",toupper(x$model$distribution))
    cat(paste0("\n",x$measure," Profile Summary : ",model))
    cat("\n\n")
    df <- as.data.frame(x$table)
    r_names <- df[,1]
    df <- df[,-1]
    rownames(df) <- r_names
    print(df, digits = digits)
    invisible(x)
}

