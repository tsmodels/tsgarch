#' Inverse AR/MA Roots
#'
#' @description Computes the inverse roots of the AR and MA polynomials of the
#' mean equation in a fitted or specified ARMA-GARCH model.
#' @param object an object of class \dQuote{tsgarch.estimate} or
#' \dQuote{tsgarch.spec}.
#' @return A list with elements \sQuote{ar} and \sQuote{ma}, each a complex
#' vector of inverse roots. The AR polynomial is
#' \eqn{1 - \phi_1 z - \dots - \phi_p z^p}{1 - phi_1 z - ... - phi_p z^p} and
#' the MA polynomial is
#' \eqn{1 + \theta_1 z + \dots + \theta_q z^q}{1 + theta_1 z + ... + theta_q z^q}.
#' If either component is absent, the corresponding vector is \code{complex(0)}.
#' @export
arma_inverse_roots <- function(object)
{
    ac <- arma_coefficients(object)
    out <- list(ar = complex(0), ma = complex(0))
    if (length(ac$ar) > 0) {
        arpoly <- c(1, -as.numeric(ac$ar))
        r <- polyroot(arpoly)
        out$ar <- 1/r
        out$ar[Mod(r) == 0] <- NaN + 1i * NaN
    }
    if (length(ac$ma) > 0) {
        mapoly <- c(1, as.numeric(ac$ma))
        r <- polyroot(mapoly)
        out$ma <- 1/r
        out$ma[Mod(r) == 0] <- NaN + 1i * NaN
    }
    return(out)
}

#' ARMA Impulse Response Function
#'
#' @description Computes the impulse response function (IRF) for the mean
#' equation of a fitted or specified ARMA-GARCH model.
#' @param object an object of class \dQuote{tsgarch.estimate} or
#' \dQuote{tsgarch.spec}.
#' @param lag.max integer, the maximum lag. If \code{NULL}, an adaptive
#' horizon is chosen by truncating when the response falls below
#' \code{1e-3} of the peak absolute value, with a floor of 10 and a cap of 100.
#' @return A list with elements \sQuote{psi} (the IRF, indexed from lag 0) and
#' \sQuote{cumulative} (the cumulative sum).
#' @export
arma_irf <- function(object, lag.max = NULL)
{
    ac <- arma_coefficients(object)
    ar <- as.numeric(ac$ar)
    ma <- as.numeric(ac$ma)
    L <- lag.max
    if (is.null(L)) {
        L <- 200
    }
    psi <- c(1, as.numeric(stats::ARMAtoMA(ar = ar, ma = ma, lag.max = L)))
    if (is.null(lag.max)) {
        threshold <- 1e-3 * max(abs(psi))
        last <- max(which(abs(psi) >= threshold))
        last <- min(max(last, 10), 100)
        psi <- psi[1:(last + 1)]
    }
    cum <- cumsum(psi)
    names(psi) <- paste0("lag", seq_along(psi) - 1)
    return(list(psi = psi, cumulative = cum))
}

#' Near AR/MA Root Cancellation
#'
#' @description Identifies AR and MA inverse roots that are close in the
#' complex plane, a visual cue for near-common-factors in the mean equation.
#' @param roots a list, normally the output of \code{\link{arma_inverse_roots}}.
#' @param tol numeric distance threshold in the complex plane.
#' @return A \code{data.table} with one row per near-cancelling pair and columns
#' \sQuote{ar_index}, \sQuote{ma_index}, \sQuote{ar_root}, \sQuote{ma_root}
#' and \sQuote{distance}.
#' @export
arma_near_cancellation <- function(roots, tol = 0.1)
{
    ar <- as.complex(roots$ar)
    ma <- as.complex(roots$ma)
    out <- data.table::data.table(ar_index = integer(0), ma_index = integer(0),
                                  ar_root = complex(0), ma_root = complex(0),
                                  distance = numeric(0))
    if (length(ar) == 0 || length(ma) == 0) return(out)
    dm <- sapply(ar, function(a) Mod(a - ma))
    if (is.matrix(dm)) {
        idx <- which(dm < tol, arr.ind = TRUE)
        if (nrow(idx) > 0) {
            for (k in seq_len(nrow(idx))) {
                i <- idx[k, 1]; j <- idx[k, 2]
                out <- rbind(out, data.table::data.table(
                    ar_index = i, ma_index = j,
                    ar_root = ar[i], ma_root = ma[j],
                    distance = dm[i, j]), use.names = TRUE)
            }
        }
    } else if (length(dm) > 0) {
        # length(ar)==1 or length(ma)==1 path, dm is a vector
        idx <- which(dm < tol)
        for (k in idx) {
            i <- ifelse(length(ar) > 1, k, 1)
            j <- ifelse(length(ma) > 1, k, 1)
            out <- rbind(out, data.table::data.table(
                ar_index = i, ma_index = j,
                ar_root = ar[i], ma_root = ma[j],
                distance = dm[k]), use.names = TRUE)
        }
    }
    return(out)
}
