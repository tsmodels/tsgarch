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
                    # NULL fallbacks cover spec/estimate objects serialized
                    # by package versions predating the xreg slot
                    "xreg" = if (!is.null(x$xreg) && !is.null(x$xreg$xreg)) xts(x$xreg$xreg, x$target$index) else xts(matrix(0, ncol = 1, nrow = length(x$target$index)), x$target$index),
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
                    "tau" = parmatrix[group == "tau"]$value,
                    "distribution" = parmatrix[group == "distribution"]$value,
                    "arpacf" = parmatrix[group == "arpacf"]$value,
                    "mapacf" = parmatrix[group == "mapacf"]$value,
                    # see arma_coefficients(): value is the target ar/ma
                    # coefficients directly when the whole polynomial is
                    # fixed, else the raw pacf parameterization
                    "ar" = if (identical(arma_block_status(parmatrix, "arpacf"), "fixed")) parmatrix[group == "arpacf"]$value else pacf_to_ar(parmatrix[group == "arpacf"]$value),
                    "ma" = if (identical(arma_block_status(parmatrix, "mapacf"), "fixed")) parmatrix[group == "mapacf"]$value else pacf_to_ma(parmatrix[group == "mapacf"]$value))

    return(value)
}

#' ARMA mean equation coefficients
#'
#' @description Returns the AR and MA coefficients of the estimated (or
#' specified) ARMA mean equation (see \code{\link{garch_modelspec}}).
#' Ordinarily these are recovered from the raw Durbin-Levinson (pacf-space)
#' parameters used internally during estimation, with stationarity of the AR
#' polynomial and invertibility of the MA polynomial guaranteed by
#' construction. If, however, an entire AR and/or MA polynomial has been
#' fixed (every \sQuote{arpacf}/\sQuote{mapacf} row has \sQuote{estimate ==
#' 0}; see \code{\link{garch_modelspec}}), the corresponding \sQuote{value}
#' entries are themselves the target ar/ma coefficients and are returned
#' as-is.
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
    ar <- numeric(0)
    ma <- numeric(0)
    if (order[1] > 0) {
        ar_raw <- parmatrix[group == "arpacf"]$value
        ar <- if (identical(arma_block_status(parmatrix, "arpacf"), "fixed")) ar_raw else pacf_to_ar(ar_raw)
        names(ar) <- paste0("ar", seq_len(order[1]))
    }
    if (order[2] > 0) {
        ma_raw <- parmatrix[group == "mapacf"]$value
        ma <- if (identical(arma_block_status(parmatrix, "mapacf"), "fixed")) ma_raw else pacf_to_ma(ma_raw)
        names(ma) <- paste0("ma", seq_len(order[2]))
    }
    list(ar = ar, ma = ma)
}
