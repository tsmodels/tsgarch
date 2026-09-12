.plot_tsgarch_estimate_arma <- function(x, which = 1:4, cumulative = FALSE,
                                         envelope = c("bartlett", "simulate"), B = 500, ...)
{
    arma_order <- if (!is.null(x$spec$model$arma)) x$spec$model$arma else c(0,0)
    if (sum(arma_order) == 0) {
        stop("\ntype = 'arma' requires a model with an ARMA mean equation (arma != c(0,0)).")
    }
    envelope <- match.arg(envelope)
    which <- as.integer(which)
    if (any(which < 1) || any(which > 4)) stop("\nwhich must be an integer vector with values in 1:4.")
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    if (length(which) == 1) {
        .plot_arma_panel(x, which[1], cumulative, envelope, B)
    } else {
        panels <- sort(unique(which))
        n <- length(panels)
        if (n == 2) {
            par(mfrow = c(1, 2))
        } else if (n >= 3) {
            par(mfrow = c(2, 2))
        }
        for (p in panels) {
            .plot_arma_panel(x, p, cumulative, envelope, B)
        }
    }
}

.plot_arma_panel <- function(x, panel, cumulative, envelope, B)
{
    switch(panel,
           "1" = .plot_arma_inverse_roots(x),
           "2" = .plot_arma_irf(x, cumulative),
           "3" = .plot_arma_acf(x, variable = "z", envelope = envelope, B = B),
           "4" = .plot_arma_acf(x, variable = "z2", envelope = envelope, B = B))
}

.plot_arma_inverse_roots <- function(x)
{
    roots <- arma_inverse_roots(x)
    near <- arma_near_cancellation(roots, tol = 0.1)
    has_canc <- nrow(near) > 0

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(2.5, 2.5, 2.5, 0.5))
    theta <- seq(0, 2 * pi, length.out = 200)
    plot(0, 0, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1,
         type = "n", xlab = "", ylab = "", main = "Inverse Roots")
    lines(cos(theta), sin(theta), col = "gray60")
    grid()
    if (length(roots$ar) > 0) {
        points(Re(roots$ar), Im(roots$ar), pch = 19, col = "steelblue")
    }
    if (length(roots$ma) > 0) {
        points(Re(roots$ma), Im(roots$ma), pch = 4, col = "coral")
    }
    legend("topright", legend = c("AR roots", "MA roots"),
           col = c("steelblue", "coral"), pch = c(19, 4), bg = "white")
    if (has_canc) {
        for (k in seq_len(nrow(near))) {
            x1 <- Re(near$ar_root[k]); y1 <- Im(near$ar_root[k])
            x2 <- Re(near$ma_root[k]); y2 <- Im(near$ma_root[k])
            segments(x1, y1, x2, y2, col = "gray30", lty = 2)
        }
        mtext("Possible common factors", side = 3, line = 0.2, cex = 0.8)
    }
}

.plot_arma_irf <- function(x, cumulative = FALSE)
{
    irf <- arma_irf(x)
    lag <- seq_along(irf$psi) - 1
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(4, 4, 2.5, 0.5))
    plot(lag, irf$psi, type = "h", lwd = 1.2, col = "steelblue",
         xlab = "Lag", ylab = expression(psi), main = "Impulse Response")
    points(lag, irf$psi, pch = 19, cex = 0.7, col = "steelblue")
    abline(h = 0, col = "gray60")
    grid()
    if (cumulative) {
        lines(lag, irf$cumulative, col = "coral", lty = 2)
        legend("topright", legend = c("IRF", "Cumulative"),
               col = c("steelblue", "coral"), lty = c(1, 2), bg = "white")
    }
}

.plot_arma_acf <- function(x, variable = c("z", "z2"), envelope = c("bartlett", "simulate"), B = 500)
{
    variable <- match.arg(variable)
    envelope <- match.arg(envelope)
    z <- as.numeric(residuals(x, standardize = TRUE))
    z <- z[!is.na(z)]
    if (variable == "z2") {
        vals <- z^2
        main <- expression(ACF(z[t]^2))
    } else {
        vals <- z
        main <- expression(ACF(z[t]))
    }
    a <- stats::acf(vals, plot = FALSE, na.action = stats::na.pass)
    lag <- as.numeric(a$lag[-1])
    acfval <- as.numeric(a$acf[-1])
    n <- length(vals)
    ci <- 1.96 / sqrt(n)

    if (envelope == "simulate") {
        stop("\nenvelope = 'simulate' is not yet implemented. Use envelope = 'bartlett'.")
    }

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(4, 4, 2.5, 0.5))
    plot(range(c(0, lag + 0.5)), range(c(min(acfval, -ci), max(acfval, ci))),
         type = "n", xlab = "Lag", ylab = "ACF", main = main)
    abline(h = 0, col = "gray60")
    abline(h = c(-ci, ci), col = "coral", lty = 2)
    lines(lag, acfval, type = "h", lwd = 1.2, col = "steelblue")
    points(lag, acfval, pch = 19, cex = 0.7, col = "steelblue")
    grid()
}
