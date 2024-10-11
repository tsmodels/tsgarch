#' Score Method
#'
#' @param x an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return The score matrix
#' @details The function returns the scores of the negative of the log likelihood
#' at the optimal solution. The scores are from the second pass of the optimizer
#' using scaling and hence represent the scaled scores. These are then rescaled
#' before returning the matrix.
#' @method estfun tsgarch.estimate
#' @aliases estfun.tsgarch.estimate
#' @rdname estfun
#' @author Alexios Galanos
#' @export
#'
estfun.tsgarch.estimate <- function(x, ...)
{
    out <- x$scaled_scores
    N <- nrow(out)
    out <- sweep(out, 2, 1/x$parameter_scale, FUN = "*")
    colnames(out) <- x$parmatrix[estimate == 1]$parameter
    return(out)
}


#' Bread Method
#'
#' @param x an object of class \dQuote{tsgarch.estimate}.
#' @param ... not currently used.
#' @return The analytic hessian of the model.
#' @method bread tsgarch.estimate
#' @aliases bread.tsgarch.estimate
#' @rdname bread
#' @author Alexios Galanos
#' @export
#'
bread.tsgarch.estimate <- function(x, ...)
{
    D <- diag(1/x$parameter_scale, nrow = length(x$parameter_scale))
    H <- t(D) %*% x$scaled_hessian %*% D
    H <- solve(H)
    return(H)
}

meat_tsgarch <- function(x, adjust = FALSE, ...)
{
    psi <- estfun(x, ...)
    k <- NCOL(psi)
    n <- NROW(psi)
    rval <- crossprod(as.matrix(psi))
    if (adjust) rval <- n/(n - k) * rval
    rownames(rval) <- colnames(rval) <- colnames(psi)
    return(rval)
}


meatHAC_tsgarch <- function(x, prewhite = FALSE, weights = NULL,  lag = NULL,
                                     kernel = c("Bartlett", "Parzen", "Quadratic Spectral",
                                                "Truncated", "Tukey-Hanning"),
          adjust = TRUE, diagnostics = FALSE, ar.method = "ols",  ...)
{
    prewhite <- as.integer(prewhite)
    umat <- estfun(x, ...)[, , drop = FALSE]
    if (is.zoo(umat)) umat <- as.matrix(coredata(umat))
    n.orig <- n <- nrow(umat)
    k <- ncol(umat)
    if (is.null(weights)) {
        if (is.null(lag)) {
            lag <- floor(bwNeweyWest(x, order.by = NULL, weights = 1, prewhite = prewhite, ar.method = ar.method,
                                     kernel = kernel[1]))
        }
        weights <- seq(1, 0, by = -(1/(lag + 1)))
    } else {
        if (length(weights) > n) {
            warning("more weights than observations, only first n used")
            weights <- weights[1:n]
        }
    }
    index <- 1:n
    umat <- umat[index, , drop = FALSE]
    if (prewhite > 0) {
        var.fit <- try(ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method))
        if (inherits(var.fit, "try-error"))
            stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
        if (k > 1) {
            D <- solve(diag(ncol(umat)) - apply(var.fit$ar, 2:3, sum))
        } else {
            D <- as.matrix(1/(1 - sum(var.fit$ar)))
        }
        umat <- as.matrix(na.omit(var.fit$resid))
        n <- n - prewhite
    }
    utu <- 0.5 * crossprod(umat) * weights[1]
    wsum <- n * weights[1]/2
    w2sum <- n * weights[1]^2/2
    if (length(weights) > 1) {
        for (ii in 2:length(weights)) {
            utu <- utu + weights[ii] * crossprod(umat[1:(n - ii + 1), , drop = FALSE], umat[ii:n, , drop = FALSE])
            wsum <- wsum + (n - ii + 1) * weights[ii]
            w2sum <- w2sum + (n - ii + 1) * weights[ii]^2
        }
    }
    utu <- utu + t(utu)
    if (adjust) utu <- n.orig/(n.orig - k) * utu
    if (prewhite > 0) utu <- crossprod(t(D), utu) %*% t(D)
    wsum <- 2 * wsum
    w2sum <- 2 * w2sum
    bc <- n^2/(n^2 - wsum)
    df <- n^2/w2sum
    rval <- utu
    if (diagnostics) attr(rval, "diagnostics") <- list(bias.correction = bc, df = df)
    return(rval)
}


