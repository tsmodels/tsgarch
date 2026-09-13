.plot_tsgarch_estimate_arma <- function(x, which = NULL, cumulative = FALSE,
                                         envelope = c("bartlett", "simulate", "parametric"),
                                         B = 500, vcov_type = "H", ...)
{
    arma_order <- if (!is.null(x$spec$model$arma)) x$spec$model$arma else c(0,0)
    if (sum(arma_order) == 0) {
        stop("\ntype = 'arma' requires a model with an ARMA mean equation (arma != c(0,0)).")
    }
    if (is.null(which)) which <- 1:4
    envelope <- match.arg(envelope)
    which <- as.integer(which)
    if (any(which < 1) || any(which > 4)) stop("\nwhich must be an integer vector with values in 1:4.")
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    # One margin/label geometry for every panel in the call. The panel helpers
    # deliberately no longer set mar themselves: they used to, with different
    # values (panel 1 reserved a narrower left margin than panels 2-4), so the
    # plot boxes had different widths and panel 1 could not line up vertically
    # with panel 3. Owning it here also keeps the boxes aligned for whatever
    # subset of panels `which` selects. Bottom margin is small because no
    # panel draws an xlab; left margin has to clear the ylab at mgp[1].
    panel_mar <- c(2.5, 3.4, 3.6, 0.8)
    panel_mgp <- c(2.2, 0.55, 0)
    # short legend labels once the panels are sub-device sized
    compact <- length(which) > 1
    if (length(which) == 1) {
        par(mar = panel_mar, mgp = panel_mgp)
        .plot_arma_panel(x, which[1], cumulative, envelope, B, vcov_type, compact)
    } else {
        panels <- sort(unique(which))
        n <- length(panels)
        if (n == 2) {
            par(mfrow = c(1, 2), mar = panel_mar, mgp = panel_mgp)
        } else if (n >= 3) {
            par(mfrow = c(2, 2), mar = panel_mar, mgp = panel_mgp)
        }
        for (p in panels) {
            .plot_arma_panel(x, p, cumulative, envelope, B, vcov_type, compact)
        }
    }
}

# Panel legends go in the top margin rather than inside the plot box. At 2x2
# size a "topright" legend is wider than the free space in the corner and ends
# up sitting on the data it is there to explain, whereas the top margin is
# empty and spans the full panel width. Anchoring off par("usr")/par("cxy")
# keeps it a fixed number of text lines above the box on any device size,
# which an `inset` fraction of the plot region would not.
.arma_panel_legend <- function(labels, col, lty = NA, pch = 1, cex = 0.9)
{
    legend(x = par("usr")[1], y = par("usr")[4] + 1.45 * par("cxy")[2],
           legend = labels, col = col, lty = lty, pch = pch, horiz = TRUE,
           bty = "n", cex = cex, xpd = NA, x.intersp = 0.4, seg.len = 1.2)
}

# titles are drawn here, not via plot(main=), so they clear the margin legend
.arma_panel_title <- function(main)
{
    mtext(main, side = 3, line = 2.2, font = par("font.main"), cex = par("cex.main") * par("cex"))
}

.plot_arma_panel <- function(x, panel, cumulative, envelope, B, vcov_type = "H", compact = FALSE)
{
    switch(panel,
           "1" = .plot_arma_inverse_roots(x, compact),
           "2" = .plot_arma_irf(x, cumulative, compact),
           "3" = .plot_arma_acf(x, variable = "z", envelope = envelope, B = B, vcov_type = vcov_type, compact = compact),
           "4" = .plot_arma_acf(x, variable = "z2", envelope = envelope, B = B, vcov_type = vcov_type, compact = compact))
}

.plot_arma_inverse_roots <- function(x, compact = FALSE)
{
    roots <- arma_inverse_roots(x)
    near <- arma_near_cancellation(roots, tol = 0.1)
    has_canc <- nrow(near) > 0

    # No par() save/restore here: this is an internal (non-exported) helper,
    # always called from .plot_tsgarch_estimate_arma(), which already
    # captures and restores the *complete* par state exactly once for the
    # whole (possibly multi-panel) call - the only place CRAN's "restore
    # graphical parameters" requirement actually applies, since this
    # function is never reachable directly by a user. A second, nested
    # full-state save/restore here would reset par("mfg") (the plotting
    # cursor within an mfrow layout) after every single panel, which is
    # exactly what previously caused every panel to overwrite the same grid
    # cell instead of advancing through the 2x2 layout.
    theta <- seq(0, 2 * pi, length.out = 200)
    # ylab is the imaginary axis, so it is worth labelling: it also fills the
    # left margin that the shared geometry reserves for panels 2-4's ylab,
    # which would otherwise be blank space here.
    plot(0, 0, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1,
         type = "n", xlab = "", ylab = "Im", main = "")
    .arma_panel_title("Inverse Roots")
    lines(cos(theta), sin(theta), col = "gray60")
    grid()
    if (has_canc) {
        for (k in seq_len(nrow(near))) {
            x1 <- Re(near$ar_root[k]); y1 <- Im(near$ar_root[k])
            x2 <- Re(near$ma_root[k]); y2 <- Im(near$ma_root[k])
            segments(x1, y1, x2, y2, col = "gray30", lty = 2)
        }
    }
    if (length(roots$ar) > 0) {
        points(Re(roots$ar), Im(roots$ar), pch = 19, col = "steelblue")
    }
    if (length(roots$ma) > 0) {
        points(Re(roots$ma), Im(roots$ma), pch = 4, col = "coral")
    }
    # Build the legend dynamically so that when only one of AR or MA is
    # present the missing component is not shown, and so that point-only
    # legend entries use lty = 0 (blank) rather than NA. legend() chokes on
    # NA line types when pch is also supplied (Error in if (do.lines) ...).
    labels <- character(0)
    col <- character(0)
    pch <- integer(0)
    lty <- integer(0)
    if (length(roots$ar) > 0) {
        labels <- c(labels, "AR roots")
        col <- c(col, "steelblue")
        pch <- c(pch, 19)
        lty <- c(lty, 0)
    }
    if (length(roots$ma) > 0) {
        labels <- c(labels, "MA roots")
        col <- c(col, "coral")
        pch <- c(pch, 4)
        lty <- c(lty, 0)
    }
    if (has_canc) {
        # the dashed joining segment explains itself as a legend key, which is
        # both clearer and cheaper on space than the mtext subtitle this
        # replaces (that subtitle also collided with the margin legend).
        labels <- c(labels, if (compact) "Common factors" else "Possible common factors")
        col <- c(col, "gray30")
        pch <- c(pch, NA)
        lty <- c(lty, 2)
    }
    .arma_panel_legend(labels, col = col, lty = lty, pch = pch)
}

.plot_arma_irf <- function(x, cumulative = FALSE, compact = FALSE)
{
    irf <- arma_irf(x)
    # psi_0 is identically 1 by the normalization of the MA(Inf)
    # representation (see arma_irf()), so it is uninformative about the fitted
    # model and only compresses the vertical scale of the lags that matter.
    # Drop it from the display, as the ACF panels do with their lag-0 spike.
    lag <- seq_along(irf$psi)[-1] - 1
    psi <- irf$psi[-1]
    # accumulate from lag 1 rather than subsetting irf$cumulative: carrying the
    # psi_0 = 1 offset would pin the line near 1 and, since it has to be inside
    # ylim to avoid being clipped, would re-compress the lags we just zoomed in
    # on. So this is the cumulative response *excluding* the contemporaneous
    # unit impact, which is on the same scale as the plotted psi.
    cum <- cumsum(psi)
    # see the comment in .plot_arma_inverse_roots(): no par() save/restore here
    # by design, and no mar either - .plot_tsgarch_estimate_arma() owns both.
    ylim <- if (cumulative) range(c(psi, cum, 0)) else range(c(psi, 0))
    plot(x = lag, psi, type = "n", ylim = ylim,
         xlab = "", ylab = expression(psi), main = "")
    .arma_panel_title("Impulse Response")
    abline(h = 0, col = "gray60")
    grid()
    lines(lag, psi, type = "h", lwd = 1.2, col = "steelblue")
    points(lag, psi, pch = 19, cex = 0.7, col = "steelblue")
    if (cumulative) {
        lines(lag, cum, col = "coral", lty = 2)
        .arma_panel_legend(c("IRF", if (compact) "Cumulative" else "Cumulative (from lag 1)"),
                           col = c("steelblue", "coral"), lty = c(1, 2),
                           pch = c(19, NA))
    }
}

.plot_arma_acf <- function(x, variable = c("z", "z2"), envelope = c("bartlett", "simulate", "parametric"),
                           B = 500, vcov_type = "H", compact = FALSE)
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

    env_lo <- env_hi <- NULL
    if (envelope %in% c("simulate", "parametric")) {
        L <- max(lag)
        if (envelope == "simulate") {
            dpars <- extract_model_values(x, object_type = "estimate", "distribution")
            distribution <- x$spec$distribution
            zsim <- matrix(rdist(distribution, n * B, mu = 0, sigma = 1,
                                 skew = dpars[1], shape = dpars[2], lambda = dpars[3]),
                           nrow = n, ncol = B)
        } else {
            zsim <- .parametric_standardized_residual_draws(x, B = B, vcov_type = vcov_type)
        }
        simvals <- if (variable == "z2") zsim^2 else zsim
        acf_sim <- apply(simvals, 2, function(col) {
            col <- col[!is.na(col)]
            as.numeric(stats::acf(col, plot = FALSE, lag.max = L, na.action = stats::na.pass)$acf[-1])
        })
        env_lo <- apply(acf_sim, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
        env_hi <- apply(acf_sim, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
    }

    # see the comment in .plot_arma_inverse_roots(): no par() save/restore here
    # by design, and no mar either - .plot_tsgarch_estimate_arma() owns both.
    ylim <- range(c(min(acfval, -ci), max(acfval, ci), env_lo, env_hi))
    plot(x = range(c(0, lag + 0.5)), y = ylim, type = "n", xlab = "", ylab = "ACF", main = "")
    .arma_panel_title(main)
    abline(h = 0, col = "gray60")
    grid()
    if (!is.null(env_lo)) {
        # "simulate" is a null-rejection band (centered near zero, like Bartlett);
        # "parametric" is a confidence band for the true ACF itself, centered on
        # (and shifted with) the point estimate, not on zero - see plot.tsgarch.estimate
        # Details. Drawn with a visually distinct style so the two are not
        # mistaken for the same kind of band.
        env_col <- if (envelope == "simulate") "darkgreen" else "purple"
        env_lty <- if (envelope == "simulate") 3 else 4
        lines(lag, env_lo, col = env_col, lty = env_lty, lwd = 1.5)
        lines(lag, env_hi, col = env_col, lty = env_lty, lwd = 1.5)
    }
    abline(h = c(-ci, ci), col = "coral", lty = 2)
    lines(lag, acfval, type = "h", lwd = 1.2, col = "steelblue")
    points(lag, acfval, pch = 19, cex = 0.7, col = "steelblue")
    # The full descriptive labels do not fit a quarter-device panel at any
    # readable cex, so `compact` trades them for short keys; the null-band vs
    # confidence-band distinction they spell out is in ?plot.tsgarch.estimate
    # and is still shown in full when a single panel owns the device.
    legend_labels <- c(if (compact) "Bartlett" else "Bartlett (null band)")
    legend_col <- c("coral")
    legend_lty <- c(2)
    if (!is.null(env_lo)) {
        if (compact) {
            env_label <- if (envelope == "simulate") "Simulated" else "Parametric CI"
        } else {
            env_label <- if (envelope == "simulate") "Simulated (null band, no param. uncertainty)" else "Parametric 95% CI of true ACF (not a null band)"
        }
        legend_labels <- c(legend_labels, env_label)
        legend_col <- c(legend_col, env_col)
        legend_lty <- c(legend_lty, env_lty)
    }
    .arma_panel_legend(legend_labels, col = legend_col, lty = legend_lty,
                       cex = if (compact) 0.9 else 0.8)
}
