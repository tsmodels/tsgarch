extract_model_values <- function(object, object_type, value_name, ...)
{
    group <- NULL
    parmatrix <- object$parmatrix
    if (object_type == "estimate") {
        x <- object$spec
    } else {
        x <- object
    }
    value <- switch(value_name,
                    "y" = x$target$y,
                    "vreg" = xts(x$vreg$vreg, x$target$index),
                    "mu" = parmatrix[group == "mu"]$value,
                    "omega" = parmatrix[group == "omega"]$value,
                    "phi" = parmatrix[group == "phi"]$value,
                    "rho" = parmatrix[group == "rho"]$value,
                    "alpha" = parmatrix[group == "alpha"]$value,
                    "gamma" = parmatrix[group == "gamma"]$value,
                    "eta" = parmatrix[group == "eta"]$value,
                    "delta" = parmatrix[group == "delta"]$value,
                    "beta" = parmatrix[group == "beta"]$value,
                    "xi" = parmatrix[group == "xi"]$value,
                    "distribution" = parmatrix[group == "distribution"]$value,
                    "arpacf" = parmatrix[group == "arpacf"]$value,
                    "mapacf" = parmatrix[group == "mapacf"]$value,
                    "ar" = pacf_to_ar(parmatrix[group == "arpacf"]$value),
                    "ma" = pacf_to_ma(parmatrix[group == "mapacf"]$value))

    return(value)
}

#' ARMA mean equation coefficients
#'
#' @description Returns the AR and MA coefficients implied by the raw
#' Durbin-Levinson (pacf-space) parameters used internally during estimation
#' (see \code{\link{garch_modelspec}}). Stationarity of the AR polynomial and
#' invertibility of the MA polynomial are guaranteed by construction and do
#' not need to be checked.
#' @param object an object of class \dQuote{tsgarch.spec} or
#' \dQuote{tsgarch.estimate}.
#' @return a list with elements \sQuote{ar} and \sQuote{ma}, each a numeric
#' vector (possibly of length zero if that component is not present).
#' @export
#'
arma_coefficients <- function(object)
{
    group <- NULL
    if (!inherits(object, c("tsgarch.spec","tsgarch.estimate"))) {
        stop("object must be of class tsgarch.spec or tsgarch.estimate")
    }
    parmatrix <- object$parmatrix
    order <- if (!is.null(object$spec)) object$spec$model$arma else object$model$arma
    if (is.null(order)) order <- c(0,0)
    ar <- if (order[1] > 0) pacf_to_ar(parmatrix[group == "arpacf"]$value) else numeric(0)
    ma <- if (order[2] > 0) pacf_to_ma(parmatrix[group == "mapacf"]$value) else numeric(0)
    if (order[1] > 0) names(ar) <- paste0("ar", seq_len(order[1]))
    if (order[2] > 0) names(ma) <- paste0("ma", seq_len(order[2]))
    list(ar = ar, ma = ma)
}
