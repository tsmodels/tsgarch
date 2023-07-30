.simulate_garch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
    }
    if (is.null(var_init)) {
        initv <- omega/(1 - sum(alpha) - sum(beta))
    } else {
        initv <- var_init
    }

    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_sqr_sim[,seq_len(maxpq)] <- matrix(initv, ncol = maxpq, nrow = nsim, byrow = TRUE)
        sigma_sim[,seq_len(maxpq)] <- matrix(sqrt(initv), ncol = maxpq, nrow = nsim, byrow = TRUE)
        epsilon[,seq_len(maxpq)] <- matrix(z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)], ncol = maxpq, nrow = nsim, byrow = TRUE)
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + alpha[j] * epsilon[,i - j]^2
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + beta[j] * sigma_sqr_sim[,i - j]
            }
        }
        sigma_sim[,i] <- sqrt(sigma_sqr_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

# adjust maxpq <- to matrix
.simulate_egarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    kappa <- egarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta)
        initv <- omega/(1 - p)
    } else {
        initv <- log(var_init)
    }
    sigma_sim <- sigma_log_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_log_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- exp(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_log_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_log_sim[,i] <- sigma_log_sim[,i] + alpha[j] * z[,i - j] + gamma[j] * (abs(z[,i - j]) - kappa)
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_log_sim[,i] <- sigma_log_sim[,i] + beta[j] * sigma_log_sim[,i - j]
            }
        }
        sigma_sim[,i] <- sqrt(exp(sigma_log_sim[,i]))
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_aparch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    if (is.null(var_init)) {
        kappa <- aparch_moment_v(distribution = distribution, gamma = gamma, delta = delta,
                                           skew = dist[1], shape = dist[2], lambda = dist[3])
        p <- sum(beta) + sum(alpha * kappa)
        initv <- omega/(1 - p)
    } else {
        initv <- var_init^(delta/2)
    }

    sigma_sim <- sigma_power_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_power_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- initv^(1/delta)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_power_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_power_sim[,i] <- sigma_power_sim[,i] + alpha[j] * (abs(epsilon[,i - j]) - gamma[j] * epsilon[,i - j])^delta
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_power_sim[,i] <- sigma_power_sim[,i] + beta[j] * sigma_power_sim[,i - j]
            }
        }
        sigma_sim[,i] <- sigma_power_sim[,i]^(1/delta)
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_gjrgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    if (is.null(var_init)) {
        kappa <- gjrgarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
        p <- sum(beta) + sum(alpha) + sum(gamma * kappa)
        initv <- omega/(1 - p)
    } else {
        initv <- var_init
    }

    sigma_sim <- sigma_squared_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_squared_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_squared_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_squared_sim[,i] <- sigma_squared_sim[,i] + alpha[j] * epsilon[,i - j]^2 + gamma[j] * (epsilon[,i - j]^2 * (epsilon[,i - j] <= 0))
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_squared_sim[,i] <- sigma_squared_sim[,i] + beta[j] * sigma_squared_sim[i - j]
            }
        }
        sigma_sim[,i] <- sqrt(sigma_squared_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_fgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    eta <- object$parmatrix[group == "eta"]$value
    beta <- object$parmatrix[group == "beta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    if (is.null(var_init)) {
        kappa <- fgarch_moment_v(distribution = distribution, gamma = gamma, eta = eta, delta = delta,
                                 skew = dist[1], shape = dist[2], lambda = dist[3])
        p <- sum(beta) + sum(alpha * kappa)
        initv <- omega/(1 - p)
    } else {
        initv <- var_init^(delta/2)
    }

    sigma_sim <- sigma_power_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_power_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- initv^(1/delta)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_power_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_power_sim[,i] <- sigma_power_sim[,i] + alpha[j] * sigma_power_sim[,i - j] * (abs(z[,i - j] - eta[j]) - gamma[j] * (z[,i - j] - eta[j]))^delta
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_power_sim[,i] <- sigma_power_sim[,i] + beta[j] * sigma_power_sim[,i - j]
            }
        }
        sigma_sim[,i] <- sigma_power_sim[,i]^(1/delta)
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_cgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    rho <- object$parmatrix[group == "rho"]$value
    phi <- object$parmatrix[group == "phi"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    if (is.null(var_init)) {
        initv <- omega/(1 - rho)
        initq <- omega/(1 - rho)
    } else {
        if (!is.matrix(var_init)) {
            if (length(var_init) != maxpq) stop(paste0("\nvar_init must be of length ", maxpq))
            initv <- var_init
            initq <- omega/(1 - rho)
        } else {
            if (nrow(var_init) != maxpq) stop(paste0("\nvar_init must have ", maxpq, " rows"))
            initv <- var_init[,1]
            initq <- var_init[,2]
        }
    }
    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    permanent_component_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        permanent_component_sim[,seq_len(maxpq)] <- initq
        sigma_sqr_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        permanent_component_sim[,i] <- variance_intercept[i] + rho * permanent_component_sim[,i - 1] + phi * (sigma_sqr_sim[,i - 1] - permanent_component_sim[,i - 1])
        sigma_sqr_sim[,i] <- permanent_component_sim[,i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + alpha[j] * (epsilon[,i - j]^2 - permanent_component_sim[,i - j])
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + beta[j] * (sigma_sqr_sim[,i - j] - permanent_component_sim[,i - j])
            }
        }
        sigma_sim[,i] <- sqrt(sigma_sqr_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    # return
    transitory_component_sim <- sigma_sim - permanent_component_sim
    sigma <- sigma_sim
    permanent_component <- permanent_component_sim
    transitory_component <- transitory_component_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq)]
        permanent_component <- permanent_component[,-seq_len(maxpq)]
        transitory_component <- transitory_component[,-seq_len(maxpq)]
        series <- series[,-seq_len(maxpq)]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    class(permanent_component) <- "tsmodel.distribution"
    class(transitory_component) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    attr(permanent_component, "date_class") <- "numeric"
    attr(transitory_component, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series, transitory_component = transitory_component, permanent_component = permanent_component)
    return(out)
}


.simulate_igarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    if (!is.null(vreg)) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        variance_intercept <- omega + vreg
        if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(0, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(0, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    if (is.null(var_init)) {
        # set the initial to some value close to integrated
        initv <- omega/0.999
    } else {
        initv <- var_init
    }

    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_sqr_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(initv)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }
    order <- object$model$order
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr_sim[,i] <- variance_intercept[i]
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + alpha[j] * epsilon[,i - j]^2
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                sigma_sqr_sim[,i] <- sigma_sqr_sim[,i] + beta[j] * sigma_sqr_sim[,i - j]
            }
        }
        sigma_sim[,i] <- sqrt(sigma_sqr_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    sigma <- sigma_sim
    series <- series_sim
    if (maxpq > 0) {
        sigma <- sigma[,-seq_len(maxpq), drop = FALSE]
        series <- series[,-seq_len(maxpq), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}
