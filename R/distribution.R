distribution_parameters <- function(distribution)
{
    tmp <- NULL
    if (distribution == "norm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "ged") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 2, lower = 0.1, upper = 50, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "std") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 100, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "snorm") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.5, lower = 0.1, upper = 10, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "sged") {
        tmp <- rbind(
            data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
            data.table(parameter = "shape", value = 2, lower = 0.1, upper = 60, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
            data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "sstd") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 1, lower = 0.01, upper = 30, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 4, lower = 2.01, upper = 60, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "nig") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 0.4, lower = 0.01, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "gh") {
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.2, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 2, lower = 0.25, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "jsu") {
        # johnson has 2 shape parameters. The second one we model with the "skew"
        # representation in rugarch
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0, lower = -20, upper = 20, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 1, lower = 0.1, upper = 10, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
    if (distribution == "ghst") {
        # setting lower bound on shape for the existence of the fourth moment
        tmp <- rbind(tmp,
                     data.table(parameter = "skew", value = 0.1, lower = -80, upper = 80, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                     data.table(parameter = "shape", value = 8.5, lower = 8.02, upper = 25, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"),
                     data.table(parameter = "lambda", value = -0.5, lower = -6, upper = 6, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\lambda"))
        return(tmp)
    }
}

distribution_class <- function(distribution)
{
    switch(distribution,
           "norm" = 1,
           "std" = 2,
           "snorm" = 3,
           "sstd" = 4,
           "ged" = 5,
           "sged" = 6,
           "nig" = 7,
           "gh" = 8,
           "jsu" = 9,
           "ghst" = 10
           )
}

distribution_abb <- function(distribution)
{
    switch(distribution,
           "norm" = "N",
           "std" = "T",
           "snorm" = "SN",
           "sstd" = "ST",
           "ged" = "GED",
           "sged" = "SGED",
           "nig" = "NIG",
           "gh" = "GH",
           "jsu" = "JSU",
           "ghst" = "GHST"
    )
}
# validation function for currently implemented distriibutions
valid_distributions <- function()
{
    c("norm", "std", "snorm", "sstd", "ged", "sged", "nig", "gh", "jsu", "ghst")
}

# Absolute Moments (egarch model)
egarch_moment <- function(distribution = "norm", skew = 0.9, shape = 4, lambda = -0.5)
{
    moment <- switch(distribution[1],
           "norm" = .norm_egarch_moment(),
           "snorm" = .snorm_egarch_moment(skew),
           "std" = .std_egarch_moment(shape),
           "sstd" = .sstd_egarch_moment(skew, shape),
           "ged"  = .ged_egarch_moment(shape),
           "sged" = .sged_egarch_moment(skew, shape),
           "jsu" = .jsu_egarch_moment(skew, shape),
           "nig" = .nig_egarch_moment(skew, shape),
           "gh" = .gh_egarch_moment(skew, shape, lambda),
           "ghst"  = .ghst_egarch_moment(skew, shape)
           )
    return(unname(moment))
}

.norm_egarch_moment <- function()
{
    moment <- sqrt(2/pi)
    return(moment)
}

.snorm_egarch_moment <- function(skew)
{
    f <- function(x) abs(x) * ddist("snorm", x, mu = 0, sigma = 1, skew = skew, shape = 0)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.ged_egarch_moment <- function(shape)
{
    lambda <- sqrt(2^(-2/shape) * gamma(1/shape)/gamma(3/shape))
    moment <- (2^(1/shape) * lambda) * gamma(2/shape)/gamma(1/shape)
    return(moment)
}

.sged_egarch_moment <- function(skew, shape)
{
    f <- function(x) abs(x) * ddist("sged", x, mu = 0, sigma = 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.std_egarch_moment <- function(shape)
{
    moment <- 2/sqrt(pi)*gamma((shape + 1)/2)/gamma(shape/2)*sqrt(shape - 2)/(shape - 1)
    return(moment)
}

.sstd_egarch_moment <- function(skew, shape)
{
    f <- function(x) abs(x) * ddist("sstd", x, mu = 0, sigma = 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.jsu_egarch_moment <- function(skew, shape)
{
    f <- function(x) abs(x) * ddist("jsu", x, mu = 0, sigma = 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.nig_egarch_moment <- function(skew, shape)
{
    f <- function(x) abs(x) * ddist("nig", x, mu = 0, sigma = 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.gh_egarch_moment <- function(skew, shape, lambda)
{
    f <- function(x) abs(x) * ddist("gh", x, mu = 0, sigma = 1, skew = skew, shape = shape, lambda = lambda)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.ghst_egarch_moment <- function(skew, shape)
{
    f <- function(x) abs(x) * ddist("ghst", x, mu = 0, sigma = 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

# Prob x<0 (gjr model)
gjrgarch_moment <- function(distribution = "norm", skew = 0.9, shape = 4, lambda = -0.5)
{
    moment <- switch(distribution[1],
                      "norm" = .norm_gjrgarch_moment(),
                      "snorm" = .snorm_gjrgarch_moment(skew),
                      "std" = .std_gjrgarch_moment(shape),
                      "sstd" = .sstd_gjrgarch_moment(skew, shape),
                      "ged"  = .ged_gjrgarch_moment(shape),
                      "sged" = .sged_gjrgarch_moment(skew, shape),
                      "jsu" = .jsu_gjrgarch_moment(skew, shape),
                      "nig" = .nig_gjrgarch_moment(skew, shape),
                      "gh" = .gh_gjrgarch_moment(skew, shape, lambda),
                      "ghst"  = .ghst_gjrgarch_moment(skew, shape)
    )
    return(unname(moment))
}


.norm_gjrgarch_moment <- function()
{
    moment <- 0.5
    return(moment)
}

.snorm_gjrgarch_moment <- function(skew)
{
    f <- function(x) ddist("snorm", x, 0, 1, skew = skew, shape = 0)
    neg_prob <- integrate(f, -Inf, 0)$value
    return(neg_prob)
}

.std_gjrgarch_moment <- function(shape)
{
    moment <- 0.5
    return(moment)
}

.sstd_gjrgarch_moment <- function(skew, shape)
{
    f <- function(x) ddist("sstd", x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

.ged_gjrgarch_moment <- function(shape)
{
    moment <- 0.5
    return(moment)
}

.sged_gjrgarch_moment <- function(skew, shape)
{
    f <- function(x) ddist("sged", x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

.jsu_gjrgarch_moment <- function(skew, shape)
{
    f <- function(x) ddist("jsu", x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

.nig_gjrgarch_moment <- function(skew, shape)
{
    f <- function(x) ddist("nig", x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

.gh_gjrgarch_moment <- function(skew, shape, lambda)
{
    f <- function(x) ddist("gh", x, 0, 1, skew = skew, shape = shape, lambda = lambda)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

.ghst_gjrgarch_moment <- function(skew, shape)
{
    f <- function(x) ddist("ghst", x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, 0)$value
    return(moment)
}

# Box-Cox transformed moment (aparch model)
aparch_moment <- function(distribution = "norm", gamma = 0.01, delta = 2, skew = 0.9, shape = 4, lambda = -0.5)
{
    moment <- switch(distribution[1],
                       "norm" = .norm_aparch_moment(gamma, delta),
                       "snorm" = .snorm_aparch_moment(gamma, delta, skew),
                       "std" = .std_aparch_moment(gamma, delta, shape),
                       "sstd" = .sstd_aparch_moment(gamma, delta, skew, shape),
                       "ged"  = .ged_aparch_moment(gamma, delta, shape),
                       "sged" = .sged_aparch_moment(gamma, delta, skew, shape),
                       "jsu" = .jsu_aparch_moment(gamma, delta, skew, shape),
                       "nig" = .nig_aparch_moment(gamma, delta, skew, shape),
                       "gh" = .gh_aparch_moment(gamma, delta, skew, shape, lambda),
                       "ghst"  = .ghst_aparch_moment(gamma, delta, skew, shape)
    )
    return(unname(moment))
}

aparch_moment_v <- Vectorize(aparch_moment, USE.NAMES = FALSE)

.norm_aparch_moment <- function(gamma, delta)
{
    moment <- 1/sqrt(pi) * ( (1 - gamma)^delta + (1 + gamma)^delta ) * 2^(0.5 * delta - 1) * gamma((delta + 1)/2)
    return(moment)
}

.snorm_aparch_moment <- function(gamma, delta, skew)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("snorm",x, 0, 1, skew = skew)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.ged_aparch_moment <- function(gamma, delta, shape)
{
    moment <- 0.5 * ((1 - gamma)^delta + (1 + gamma)^delta) * gamma(1/shape)^(0.5 * delta - 1) * gamma(3/shape)^(-0.5 * delta) * gamma((1 + delta)/shape)
    return(moment)
}

.sged_aparch_moment <- function(gamma, delta, skew, shape)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("sged",x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.std_aparch_moment <- function(gamma, delta, shape)
{
    moment <- ((shape - 2)^(delta/2) * gamma((shape - delta)/2) * gamma((delta + 1)/2))/(gamma(shape/2) * 2 * sqrt(pi)) * ((1 - gamma)^delta + (1 + gamma)^delta)
    return(moment)
}

.sstd_aparch_moment <- function(gamma, delta, skew, shape)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("sstd",x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.jsu_aparch_moment <- function(gamma, delta, skew, shape)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("jsu",x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.nig_aparch_moment <- function(gamma, delta, skew, shape)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("nig",x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.gh_aparch_moment <- function(gamma, delta, skew, shape, lambda)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("gh",x, 0, 1, skew = skew, shape = shape, lambda = lambda)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.ghst_aparch_moment <- function(gamma, delta, skew, shape)
{
    f <- function(x) (abs(x) - gamma*x)^delta * ddist("ghst",x, 0, 1, skew = skew, shape = shape)
    moment <- integrate(f, -Inf, Inf)$value
    return(moment)
}

.norm_aparch_moment_jacobian <- function(alpha, gamma, delta)
{
    alpha_grad <- (2^(delta/2 - 1) * ((1 - gamma)^delta + (1 + gamma)^delta) * gamma((1 + delta)/2))/sqrt(pi)
    gamma_grad <- (2^(delta/2 - 1) * alpha * (-(1 - gamma)^(delta - 1) + (1 + gamma)^(delta - 1)) * delta * gamma((1 + delta)/2))/sqrt(pi)
    delta_grad <- (2^(delta/2 - 2) * alpha * (gamma((delta + 1)/2) * ((1 - gamma)^delta * log(2) + (1 + gamma)^delta * log(2)  + 2 * (1 - gamma)^delta * log(1 - gamma) +
                                                                          2 * (1 + gamma)^delta * log(gamma + 1) + ((1 - gamma)^delta + (1 + gamma)^delta) * psigamma((delta + 1)/2, 0))))/sqrt(pi)
    out <- c(alpha_grad, gamma_grad, delta_grad)
    return(out)
}

.std_aparch_moment_jacobian <- function(alpha, gamma, delta, shape)
{
    denonimator1 <- (2 * sqrt(pi) * gamma(shape/2))
    denonimator2 <- (4 * sqrt(pi) * gamma(shape/2))
    de1 <- gamma(1/2 * (shape - delta))
    de2 <- gamma((1 + delta)/2)
    de3 <- ((1 - gamma)^delta + (1 + gamma)^delta)
    dp <- delta/2
    alpha_grad <- (de3 * (shape - 2)^(delta/2) * gamma((1 + delta)/2) * gamma(1/2 * (-delta + shape))) / denonimator1
    gamma_grad <- (alpha * delta * ((1 + gamma)^(delta - 1) - (1 - gamma)^(delta - 1)) * (shape - 2)^dp * de2 * de1) / denonimator1
    delta_grad <- (alpha/denonimator2) * (shape - 2)^dp * de2 * de1 *
        (((1 - gamma)^delta * (2 * log(1 - gamma) + log(shape - 2))) + ((1 + gamma)^delta * (2 * log(1 + gamma) + log(shape - 2)))  +
        (de3 * (psigamma((1 + delta)/2, 0) - psigamma((shape - delta)/2))))
    shape_grad <- (alpha/denonimator2) * de3 * (shape - 2)^(delta/2 - 1) * de2 * de1 * (delta - (shape - 2) * psigamma(shape/2, 0) + (shape - 2) * psigamma((shape - delta)/2, 0))
    out <- c(alpha_grad, gamma_grad, delta_grad, shape_grad)
    return(out)
}

.ged_aparch_moment_jacobian <- function(alpha, gamma, delta, shape)
{
    de1 <- gamma(1/shape)^(0.5 * delta - 1) * gamma(3/shape)^(-0.5 * delta) * gamma((1 + delta)/shape)
    alpha_grad <- 0.5 * ((1 - gamma)^delta + (1 + gamma)^delta) * de1
    gamma_grad <- 0.5 * alpha * (-1 * (1 - gamma)^(delta - 1) + (1 + gamma)^(delta - 1)) * delta * de1
    delta_grad <- (1/shape) * alpha * de1 *
        (shape * (0.5 * (1 - gamma)^delta * log(1 - gamma) + 0.5 * (1 + gamma)^delta * log(1 + gamma) + ((1 - gamma)^delta + (1 + gamma)^delta) * (0.25 * lgamma(1/shape) - 0.25 * lgamma(3/shape))) +
             (0.5 * (1 - gamma)^delta + 0.5 * (1 + gamma)^delta) * psigamma((1 + delta)/shape,0))
    shape_grad <- (1/(shape^2)) * alpha * (-0.25 * (1 - gamma)^delta - 0.25 * (1 + gamma)^delta) * de1 *
        ((delta - 2) * psigamma(1/shape, 0) - 3 * delta * psigamma(3/shape, 0) + (2 * delta + 2) * psigamma((1 + delta)/shape,0))
    out <- c(alpha_grad, gamma_grad, delta_grad, shape_grad)
    return(out)
}


fgarch_moment <- function(distribution = "norm", gamma = 0.01, eta = 0.0, delta = 2, skew = 0.9, shape = 4, lambda = -0.5)
{
    f <- function(x) (abs(x - eta) - gamma * (x - eta))^delta * ddist(distribution, x, mu = 0, sigma = 1, skew = skew, shape = shape, lambda = lambda)
    moment <- integrate(f, -Inf, Inf)$value
    moment <- unname(moment)
    return(moment)
}

fgarch_moment_v <- Vectorize(fgarch_moment, USE.NAMES = FALSE)
