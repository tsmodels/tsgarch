.distribution_abb <- function(distribution, standardized = FALSE)
{
    out <- switch(distribution,
               "norm" = "N\\left(0,1\\right)",
               "snorm" = "SN\\left(0,1,\\zeta\\right)",
               "std" = "T\\left(0,1,\\nu\\right)",
               "sstd" = "ST\\left(0,1,\\zeta, \\nu\\right)",
               "ged" = "GED\\left(0,1,\\nu\\right)",
               "sged" = "SGED\\left(0,1,\\zeta, \\nu\\right)",
               "jsu" = "JSU\\left(0,1,\\zeta, \\nu\\right)",
               "nig" = "NIG\\left(0,1,\\zeta, \\nu\\right)",
               "ghst" = "GHST\\left(0,1,\\zeta, \\nu\\right)",
               "ghyp" = "GH\\left(0,1,\\zeta, \\nu,\\lambda\\right)")
    if (!standardized) {
        out <- gsub("0,1","0,\\\\sigma_t",out)
    }
    return(out)
}

.equation_regressors <- function(vreg = NULL, variance_targeting = FALSE)
{
    eq2 <- eq1 <- NULL
    if (!is.null(vreg)) {
        n <- NCOL(vreg)
        if (n == 1) {
            eq1 <- "\\xi_1 v_{1,t}"
            if (variance_targeting) {
                eq2 <- "\\xi_1 \\hat v_1"
            }
        } else {
            eq1 <- paste0("\\sum_{j=1}^",n," \\xi_j v_{j,t}")
            if (variance_targeting) {
                eq2 <- paste0("\\sum_{j=1}^",n," \\xi_j \\hat v_{j}")
            }
        }
    }
    return(list(regressor = eq1, variance_target_regressor = eq2))
}


# garch equations
.equation_garch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat\\sigma^2 (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")
    eq_distribution <- .distribution_abb(distribution)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\varepsilon^2_{t-j}")
        } else {
            eq_alpha <- paste0("\\alpha_1\\varepsilon^2_{t-1}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\sigma^2_{t-j}")
        } else {
            eq_beta <- paste0("\\beta_1\\sigma^2_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    # persistence
    eq_persistence <- paste0("P = \\sum_{j=1}^q \\alpha_j + \\sum_{j=1}^p \\beta_j")
    if (multiplicative) {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\frac{\\exp(\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j})}{1 - P}")
    } else {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - P}")
    }
    # collect equations
    eq_garch <- paste0("\\sigma^2_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)
    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}

# egarch equations
.equation_egarch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "log(\\hat \\sigma^2) (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        m <- ncol(vreg)
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    eq_distribution <- .distribution_abb(distribution, standardized = TRUE)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\left(\\sum_{j=1}^",order[1],"\\alpha_j z_{t-j} + \\gamma_j \\left(\\left|z_{t-j}\\right| - E\\left[\\left|z\\right|\\right]\\right)\\right)")
        } else {
            eq_alpha <- paste0("\\alpha_1 z_{t-1} + \\gamma_1 \\left(\\left|z_{t-1}\\right| - E\\left[\\left|z\\right|\\right]\\right)")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j log\\left(\\sigma^2_{t-j}\\right)")
        } else {
            eq_beta <- paste0("\\beta_1 log\\left(\\sigma^2_{t-1}\\right)")
        }
    } else {
        eq_beta <- NULL
    }
    eq_persistence <- paste0("P = \\sum_{j=1}^p \\beta_j")
    eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = exp\\left(\\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - P}\\right)")
    # collect equations
    eq_garch <- paste0(" log\\left(\\sigma^2_t\\right) = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)

    eq_distribution <- paste0("z_t = \\frac{\\varepsilon_t}{\\sigma_t} \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}


# aparch equations
.equation_aparch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat \\sigma^{\\delta} (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        m <- ncol(vreg)
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")

    eq_distribution <- .distribution_abb(distribution, standardized = FALSE)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\left(\\left|\\varepsilon_{t-j}\\right| - \\gamma_j \\varepsilon_{t-j}\\right)^{\\delta}")
        } else {
            eq_alpha <- paste0("\\alpha_1 \\left(\\left|\\varepsilon_{t-1}\\right| - \\gamma_1 \\varepsilon_{t-1}\\right)^{\\delta}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\sigma^{\\delta}_{t-j}")
        } else {
            eq_beta <- paste0("\\beta_1 \\sigma^{\\delta}_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    eq_persistence <- paste0("P = \\sum_{j=1}^p \\beta_j + \\sum_{j=1}^q \\alpha_j \\kappa_j,\\quad \\kappa_j = E\\left[\\left(\\left|z_{t-j}\\right| - \\gamma_j z_{t-j}\\right)^{\\delta}\\right]")
    if (multiplicative) {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\exp(\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j})}{1 - P}\\right)^{2/\\delta}")
    } else {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - P}\\right)^{2/\\delta}")
    }
    # collect equations
    eq_garch <- paste0(" \\sigma^{\\delta}_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)

    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}

# aparch equations
.equation_gjrgarch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat \\sigma^2 (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        m <- ncol(vreg)
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")

    eq_distribution <- .distribution_abb(distribution, standardized = FALSE)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\varepsilon^2_{t-j} + \\gamma_j I_{t-j} \\varepsilon^2_{t-j}")
        } else {
            eq_alpha <- paste0("\\alpha_1 \\varepsilon^2_{t-1} + \\gamma_1 I_{t-1} \\varepsilon^2_{t-1}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\sigma^2_{t-j}")
        } else {
            eq_beta <- paste0("\\beta_1 \\sigma^2_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    eq_persistence <- paste0("P = \\sum_{j=1}^p \\beta_j + \\sum_{j=1}^q \\alpha_j  + \\sum_{j=1}^q \\gamma_j \\kappa,\\quad \\kappa = E\\left[I_t z^2_t\\right]")
    if (multiplicative) {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\exp(\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j})}{1 - P}\\right)")
    } else {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - P}\\right)")
    }
    # collect equations
    eq_garch <- paste0(" \\sigma^2_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)
    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}

# fgarch equations
.equation_fgarch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat \\sigma^{\\delta} (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        m <- ncol(vreg)
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")

    eq_distribution <- .distribution_abb(distribution, standardized = TRUE)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\sigma^{\\delta}_{t-j}\\left(\\left|z_{t-j} - \\eta_j\\right| - \\gamma_j \\left(z_{t-j}-\\eta_j\\right)\\right)^{\\delta}")
        } else {
            eq_alpha <- paste0("\\alpha_1 \\sigma^{\\delta}_{t-j}\\left(\\left|z_{t-1} - \\eta_1\\right| - \\gamma_1 \\left(z_{t-1} - \\eta_1\\right)\\right)^{\\delta}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\sigma^{\\delta}_{t-j}")
        } else {
            eq_beta <- paste0("\\beta_1 \\sigma^{\\delta}_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    eq_persistence <- paste0("P = \\sum_{j=1}^p \\beta_j + \\sum_{j=1}^q \\alpha_j \\kappa_j,\\quad \\kappa_j = E\\left[\\left(\\left|z_{t-j} - \\eta_j\\right| - \\gamma_j \\left(z_{t-j} - \\eta_j\\right)\\right)^{\\delta}\\right]")
    if (multiplicative) {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\exp(\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j})}{1 - P}\\right)^{2/\\delta}")
    } else {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\left(\\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - P}\\right)^{2/\\delta}")
    }
    # collect equations
    eq_garch <- paste0(" \\sigma^{\\delta}_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)
    eq_distribution <- paste0("z_t = \\frac{\\varepsilon_t}{\\sigma_t} \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}

# garch equations
.equation_cgarch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat\\sigma^2 (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")

    eq_distribution <- .distribution_abb(distribution)
    eq_permanent_component <- paste0("q_t = ",eq_constant," + \\rho q_{t-1} + \\phi \\left(\\varepsilon^2_{t-1} - \\sigma^2_{t-1}\\right)")
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\left(\\varepsilon^2_{t-j} - q_{t-j}\\right)")
        } else {
            eq_alpha <- paste0("\\alpha_1\\left(\\varepsilon^2_{t-1}-q_{t-1}\\right)")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\left(\\sigma^2_{t-j} - q_{t-j}\\right)")
        } else {
            eq_beta <- paste0("\\beta_1\\left(\\sigma^2_{t-1}-q_{t-1}\\right)")
        }
    } else {
        eq_beta <- NULL
    }
    # persistence
    eq_persistence <- paste0("P = \\sum_{j=1}^q \\alpha_j + \\sum_{j=1}^p \\beta_j")
    if (multiplicative) {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\frac{\\exp(\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j})}{1 - \\rho}")
    } else {
        eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\frac{\\omega + \\sum_{j=1}^r \\xi_j \\hat v_{j}}{1 - \\rho}")
    }
    # collect equations
    eq_garch <- paste0("\\sigma^2_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)
    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}


.equation_igarch <- function(order, vreg = NULL, multiplicative = FALSE, distribution = "norm", variance_targeting = FALSE)
{
    regressor_equations <- .equation_regressors(vreg, variance_targeting = variance_targeting)
    # constant
    if (variance_targeting) {
        eq_constant <- "\\hat\\sigma^2 (1 - P)"
        if (!is.null(vreg)) eq_constant <- paste0("\\left(",eq_constant," - ",regressor_equations$variance_target_regressor,"\\right)")
    } else {
        eq_constant <- "\\omega"
    }
    if (!is.null(vreg)) {
        eq_constant <- paste0(eq_constant," + ", regressor_equations$regressor)
    }
    if (multiplicative) eq_constant <- paste0("\\exp(",eq_constant,")")

    eq_distribution <- .distribution_abb(distribution)
    if (order[1] > 0) {
        if (order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",order[1],"\\alpha_j \\varepsilon^2_{t-j}")
        } else {
            eq_alpha <- paste0("\\alpha_1\\varepsilon^2_{t-1}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (order[2] > 0) {
        if (order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",order[2],"\\beta_j \\sigma^2_{t-j}")
        } else {
            eq_beta <- paste0("\\beta_1\\sigma^2_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    # persistence
    eq_persistence <- paste0("P = \\sum_{j=1}^q \\alpha_j + \\sum_{j=1}^p \\beta_j == 1")
    eq_unconditional <- paste0("E\\left[\\varepsilon^2_t\\right] = \\infty")
    # collect equations
    eq_garch <- paste0("\\sigma^2_t = ",eq_constant)
    if (!is.null(eq_alpha)) eq_garch <- paste0(eq_garch, " + ", eq_alpha)
    if (!is.null(eq_beta)) eq_garch <- paste0(eq_garch, " + ", eq_beta)
    eq_distribution <- paste0("\\varepsilon_t \\sim ", eq_distribution)
    out <- list(eq_distribution = eq_distribution, eq_garch = eq_garch,eq_persistence = eq_persistence,eq_unconditional = eq_unconditional)
    return(out)
}
