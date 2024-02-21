.newsimpact_garch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) omega <- exp(omega)
    alpha <- beta <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    news_impact_curve <- omega + beta*long_run_variance + alpha*(epsilon^2)
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}


.newsimpact_egarch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    alpha <- beta <- gamma <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[1] > 0) gamma <- object$parmatrix[group == "gamma"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    kappa <- egarch_moment(distribution = object$spec$distribution, skew = skew, shape = shape, lambda = lambda)
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    long_run_sigma <- sqrt(long_run_variance)
    std_residuals <- (epsilon)/long_run_sigma
    news_impact_curve <- exp(omega + alpha * std_residuals + gamma * (abs(std_residuals) - kappa) + beta * log(long_run_variance))
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}



.newsimpact_aparch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) omega <- exp(omega)
    alpha <- beta <- gamma <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[1] > 0) gamma <- object$parmatrix[group == "gamma"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    delta <- object$parmatrix[parameter == "delta"]$value[1]
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    news_impact_curve <- omega + alpha[1] * (abs(epsilon) - gamma[1] * epsilon)^(delta)  + beta[1] * (long_run_variance^(delta/2))
    news_impact_curve <- news_impact_curve^(2/delta)
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}

.newsimpact_gjrgarch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) omega <- exp(omega)

    alpha <- beta <- gamma <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[1] > 0) gamma <- object$parmatrix[group == "gamma"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    news_impact_curve <- omega + alpha[1] * epsilon^2 + gamma[1] * (epsilon^2 * (epsilon <= 0)) + beta[1] * long_run_variance
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}



.newsimpact_fgarch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) omega <- exp(omega)

    alpha <- beta <- gamma <- eta <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[1] > 0) gamma <- object$parmatrix[group == "gamma"]$value[1]
    if (order[1] > 0) eta <- object$parmatrix[group == "eta"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    delta <- object$parmatrix[parameter == "delta"]$value[1]
    skew <- object$parmatrix[parameter == "skew"]$value
    shape <- object$parmatrix[parameter == "shape"]$value
    lambda <- object$parmatrix[parameter == "lambda"]$value
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    long_run_sigma <- sqrt(long_run_variance)
    std_residuals <- (epsilon)/long_run_sigma
    news_impact_curve <- omega + alpha * long_run_variance^(delta/2) * (abs(std_residuals - eta) - gamma * (std_residuals - eta))^(delta)  + beta * (long_run_variance^(delta/2))
    news_impact_curve <- news_impact_curve^(2/delta)
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}

.newsimpact_cgarch <- function(object, epsilon = NULL) {
    parameter <- value <- group <- lower <- upper <- NULL
    if (is.null(epsilon)) {
        epsilon_minimum <- min(residuals(object))
        if (abs(epsilon_minimum) < 1) {
            epsilon <- seq(round(epsilon_minimum, 2), round(abs(epsilon_minimum), 2), length.out = 101)
        } else{
            epsilon <- seq(round(epsilon_minimum, 0), round(abs(epsilon_minimum), 0), length.out = 101)
        }
    } else {
        epsilon <- sort(epsilon)
    }
    order <- object$spec$model$order
    omega <- omega(object)
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        omega <- omega + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) omega <- exp(omega)

    alpha <- beta <- 0
    if (order[1] > 0) alpha <- object$parmatrix[group == "alpha"]$value[1]
    if (order[2] > 0) beta <- object$parmatrix[group == "beta"]$value[1]
    rho <- object$parmatrix[group == "rho"]$value
    phi <- object$parmatrix[group == "phi"]$value
    long_run_variance <- rep(as.numeric(unconditional(object)), length(epsilon))
    permanent_component <- omega + rho * long_run_variance + phi * (epsilon^2 - long_run_variance)
    news_impact_curve <- permanent_component + beta * (long_run_variance - omega) + alpha*(epsilon^2 - omega)
    yexpr <- expression(sigma[t]^2)
    xexpr <- expression(epsilon[t - 1])
    out <- list(y = as.numeric(news_impact_curve), x = epsilon, yexpr = yexpr, xexpr = xexpr, model = toupper(object$spec$model$model))
    class(out) <- "tsgarch.newsimpact"
    return(out)
}

.newsimpact_igarch <- function(object, epsilon = NULL) {
    warning("\nnewsimpact for igarch model not available (infinite unconditional variance)")
    return(NULL)
}
