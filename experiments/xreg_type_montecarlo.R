# Monte Carlo: what does choosing the wrong xreg_type cost, in BOTH directions?
#
# Data are generated directly here (not via the package's simulate()), so the
# DGP does not inherit anything from the implementation being tested.
#
#   armax truth:       y_t = mu + phi*(y_{t-1} - mu) + tau*x_t + e_t
#   arma_errors truth: y_t = mu + tau*x_t + w_t,  w_t = phi*w_{t-1} + e_t
#
# Both are restrictions of the ADL  y_t = c + phi*y_{t-1} + b0*x_t + b1*x_{t-1} + e_t
#   armax:       b1 = 0
#   arma_errors: b1 = -phi*b0
# so each convention, fitted to the other's data, omits a term collinear with
# a different parameter. The question is where the bias lands in each case.
suppressMessages(library(tsgarch)); suppressMessages(library(xts)); suppressMessages(library(data.table))

MU <- 0.5; PHI <- 0.6; TAU <- 1.5; SIG <- 1; N <- 2000

sim_x <- function(n, rho) {
    if (rho == 0) return(rnorm(n))
    as.numeric(arima.sim(n = n, list(ar = rho)))
}
gen <- function(dgp, rho, n = N) {
    x <- sim_x(n, rho)
    e <- rnorm(n, 0, SIG)
    y <- numeric(n)
    if (dgp == "armax") {
        for (i in 2:n) y[i] <- MU + PHI * (y[i - 1] - MU) + TAU * x[i] + e[i]
    } else {
        w <- numeric(n)
        for (i in 2:n) w[i] <- PHI * w[i - 1] + e[i]
        y <- MU + TAU * x + w
    }
    keep <- 501:n
    d <- as.Date(seq_along(keep), origin = "1970-01-01")
    list(y = xts(y[keep], d), x = xts(matrix(x[keep], ncol = 1), d))
}
fit1 <- function(dat, type) {
    spec <- garch_modelspec(dat$y, model = "garch", order = c(0,0), constant = TRUE,
                            arma = c(1,0), xreg = dat$x, xreg_type = type)
    m <- try(suppressWarnings(estimate(spec)), silent = TRUE)
    if (inherits(m, "try-error")) return(c(tau = NA, phi = NA, ll = NA))
    c(tau = m$parmatrix[parameter == "tau1"]$value,
      phi = as.numeric(arma_coefficients(m)$ar),
      ll  = as.numeric(logLik(m)))
}

R <- as.integer(Sys.getenv("MCREP", "100"))
set.seed(20260914)
grid <- expand.grid(dgp = c("armax","arma_errors"), rho = c(0, 0.9), stringsAsFactors = FALSE)
res <- list()
for (g in seq_len(nrow(grid))) {
    dgp <- grid$dgp[g]; rho <- grid$rho[g]
    for (r in seq_len(R)) {
        dat <- gen(dgp, rho)
        a <- fit1(dat, "arma_errors"); b <- fit1(dat, "armax")
        res[[length(res) + 1]] <- data.table(dgp = dgp, rho = rho, rep = r,
                                             fit = c("arma_errors","armax"),
                                             tau = c(a["tau"], b["tau"]),
                                             phi = c(a["phi"], b["phi"]),
                                             ll  = c(a["ll"],  b["ll"]))
    }
    cat("done:", dgp, "rho =", rho, "\n"); flush.console()
}
res <- rbindlist(res)
saveRDS(res, "/tmp/xreg_type_mc.rds")

cat(sprintf("\nreplications per cell: %d ; failures: %d\n", R, sum(is.na(res$tau))))
cat(sprintf("true values: mu = %.2f, phi = %.2f, tau = %.2f, n = %d\n\n", MU, PHI, TAU, length(gen("armax",0)$y)))

summ <- res[!is.na(tau), .(
    tau_bias = mean(tau) - TAU, tau_rmse = sqrt(mean((tau - TAU)^2)),
    phi_bias = mean(phi) - PHI, phi_rmse = sqrt(mean((phi - PHI)^2))
), by = .(dgp, rho, fit)][order(dgp, rho, fit)]
summ[, correct := (dgp == fit)]
print(summ, digits = 3)

cat("\n--- likelihood: can it tell you which convention is right? ---\n")
wide <- dcast(res[!is.na(ll)], dgp + rho + rep ~ fit, value.var = "ll")
ll <- wide[, .(mean_ll_gap_correct_minus_wrong =
                   mean(ifelse(dgp == "armax", armax - arma_errors, arma_errors - armax)),
               pct_correct_model_higher_ll =
                   100 * mean(ifelse(dgp == "armax", armax > arma_errors, arma_errors > armax))),
           by = .(dgp, rho)][order(dgp, rho)]
print(ll, digits = 4)
