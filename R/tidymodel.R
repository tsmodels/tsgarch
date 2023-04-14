#' Tidy an estimate GARCH object
#'
#' @description Tidy for tsgarch
#' @param x an object of class \dQuote{tsgarch.estimate}.
#' @param conf.int Logical indicating whether or not to include a confidence
#' interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence interval if
#' conf.int = TRUE. Must be strictly greater than 0 and less than 1. Defaults to
#' 0.95, which corresponds to a 95 percent confidence interval.
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param vcov_type the vcov method to use.
#' @param ... not currently used.
#' @return A tibble() with columns:
#' \itemize{
#' \item{conf.high Upper bound on the confidence interval for the estimate.}
#' \item{conf.low Lower bound on the confidence interval for the estimate.}
#' \item{estimate The estimated value of the regression term.}
#' \item{p.value The two-sided p-value associated with the observed statistic.}
#' \item{statistic The value of a T-statistic to use in a hypothesis that the
#' regression term is non-zero.}
#' \item{std.error The standard error of the regression term.}
#' \item{term The name of the regression term.}
#' }
#' @aliases tidy
#' @method tidy tsgarch.estimate
#' @rdname tidy
#' @export
#'
#'
tidy.tsgarch.estimate <- function(x, conf.int = FALSE, conf.level = 0.95, vcov_type = "H", ...) {
    Estimate <- `t value` <- `Std. Error` <- `Pr(>|t|)` <- NULL
    result <- summary(x, vcov_type = vcov_type)$coefficients |>
        rename(estimate = Estimate,
                      std.error = `Std. Error`,
                      statistic = `t value`,
                      p.value = `Pr(>|t|)`)
    if (conf.int) {
        ci <- as_tibble(as.data.table(confint(x, level = conf.level), keep.rownames = TRUE) |> setnames("rn","term"), rownames = NA)
        result <- left_join(result, ci, by = "term")
    }
    result
}

#' Construct a summary glance of an estimated GARCH object
#'
#' @description glance method for tsgarch
#' @param x an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return A tibble().
#' @aliases glance
#' @method glance tsgarch.estimate
#' @rdname glance
#' @export
#'
#'
glance.tsgarch.estimate <- function(x, ...) {
    # include loglik, nobs, convergence, init variance etc, AIC, BIC
    out <- tibble(
            sigma = sqrt(unconditional(x)),
            logLik = as.numeric(logLik(x)),
            AIC = AIC(x),
            BIC = BIC(x),
            nobs = nobs(x)
        )
    return(out)
}
