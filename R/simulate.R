.simulate_garch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, burn = 0, ...)
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
    extra_args <- list(...)
    h <- h + burn
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h + burn.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x (h + burn).")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }

    if (is.null(var_init)) {
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        p <- sum(alpha) + sum(beta)
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init
    }

    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
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

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        if (length(init) != maxpq) init <- rep(init[1], maxpq)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    } else {
        init <- epsilon[,seq_len(maxpq), drop = FALSE]^2
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    }
    simc <- .garchsimvec(z = z, epsilon = epsilon, sigma_sqr_sim = sigma_sqr_sim, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, beta = beta, mu = mu)
    sigma <- simc$sigma
    series <- simc$series
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_egarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, burn = 0, ...)
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
    h <- h + burn
    extra_args <- list(...)
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
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
    kappa <- egarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta)
        initv <- (omega + mean_vreg)/(1 - p)
    } else {
        initv <- log(var_init)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }

    sigma_sim <- sigma_log_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_log_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- sqrt(exp(initv))
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]
    }

    order <- object$model$order

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        if (length(init) != maxpq) init <- rep(init[1], maxpq)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    } else {
        init <- (abs(z[,seq_len(order[1]), drop = FALSE]) - kappa)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    }

    simc <- .egarchsimvec(z = z, sigma_log_sim = sigma_log_sim, variance_intercept = variance_intercept, init = init, alpha = alpha, gamma = gamma, beta = beta, kappa = kappa, mu = mu, order = order)
    sigma <- simc$sigma
    series <- simc$series
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_aparch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                            vreg = NULL, burn = 0, ...)
{
    h <- h + burn
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
    extra_args <- list(...)
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

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

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

    kappa <- aparch_moment_v(distribution = distribution, gamma = gamma, delta = delta,
                             skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
    } else {
        initv <- var_init^(delta/2)
    }
    order <- object$model$order
    sigma_sim <- sigma_power_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    series_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    epsilon <- matrix(0, nrow = nsim, ncol = maxpq + h)
    if (maxpq > 0) {
        sigma_power_sim[,seq_len(maxpq)] <- initv
        sigma_sim[,seq_len(maxpq)] <- initv^(1/delta)
        epsilon[,seq_len(maxpq)] <- z[,seq_len(maxpq)] * sigma_sim[,seq_len(maxpq)]

        if (!is.null(extra_args$arch_initial)) {
            init <- extra_args$arch_initial
            if (length(init) != maxpq) init <- rep(init[1], maxpq)
            init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
        } else {
            if (is.null(innov_init)) {
                init <- kappa * (initv^(delta/2))
                if (length(init) != maxpq) init <- rep(init[1], maxpq)
                init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
            } else {
                init <- (abs(epsilon[,seq_len(order[1])]) - gamma * epsilon[,seq_len(order[1])])^delta
                if (length(init) != maxpq) init <- rep(init[1], maxpq)
                init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
            }
        }
    }

    simc <- .aparchsimvec(epsilon = epsilon, sigma_power_sim = sigma_power_sim, z = z, variance_intercept = variance_intercept,
                          init = init, alpha = alpha, gamma = gamma, beta = beta, delta = delta, mu = mu, order = order)
    sigma <- simc$sigma
    series <- simc$series
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_gjrgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, burn = 0, ...)
{
    h <- h + burn
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    extra_args <- list(...)
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    gamma <- object$parmatrix[group == "gamma"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }
    kappa <- gjrgarch_moment(distribution = distribution, skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha) + sum(gamma * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
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

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        if (length(init) != maxpq) init <- rep(init[1], maxpq)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    } else {
        init <- (epsilon[,seq_len(maxpq), drop = FALSE]^2 * kappa)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    }

    order <- object$model$order

    simc <- .gjrsimvec(epsilon = epsilon, sigma_sqr_sim = sigma_squared_sim, z = z, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, gamma = gamma, beta = beta, mu = mu)

    sigma <- simc$sigma
    series <- simc$series

    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series, epsilon = epsilon)
    return(out)
}

.simulate_fgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, burn = 0, ...)
{
    h <- h + burn
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
    extra_args <- list(...)
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
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

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

    kappa <- fgarch_moment_v(distribution = distribution, gamma = gamma, eta = eta, delta = delta,
                             skew = dist[1], shape = dist[2], lambda = dist[3])
    if (is.null(var_init)) {
        p <- sum(beta) + sum(alpha * kappa)
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - p)
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

    if (maxpq > 0) {
        if (!is.null(extra_args$arch_initial)) {
            init <- extra_args$arch_initial
            if (length(init) != maxpq) init <- rep(init[1], maxpq)
            init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
        } else {
            if (is.null(innov_init)) {
                init <- kappa
                init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
            } else {
                init <- (abs(z[,seq_len(maxpq)] - eta) - gamma * (z[,seq_len(maxpq)] - eta))^delta
                init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
            }
        }
    }

    order <- object$model$order

    simc <- .fgarchsimvec(epsilon = epsilon, sigma_power_sim = sigma_power_sim, z = z, variance_intercept = variance_intercept, init = init,
                          alpha = alpha, gamma = gamma, eta = eta, beta = beta, delta = delta, mu = mu, order = order)
    sigma <- simc$sigma
    series <- simc$series
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}

.simulate_cgarch <- function(object, h = 1000, seed  = NULL, nsim = 1, var_init = NULL, innov = NULL, innov_init = NULL,
                             vreg = NULL, burn = 0, ...)
{
    h <- h + burn
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
    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x h.")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }
    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- innov_init
    }

    if (is.null(var_init)) {
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - rho)
        initq <- numerator/(1 - rho)
        if (object$vreg$multiplicative) {
            initv <- exp(initv)
            initq <- exp(initq)
        }
    } else {
        if (!is.matrix(var_init)) {
            if (length(var_init) != maxpq) stop(paste0("\nvar_init must be of length ", maxpq))
            initv <- var_init
            initq <- var_init
        } else {
            if (nrow(var_init) != maxpq) stop(paste0("\nvar_init must have ", maxpq, " rows"))
            initq <- var_init[,1]
            initv <- var_init[,2]
        }
    }
    sigma_sim <- sigma_sqr_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    permanent_component_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
    transitory_component_sim <- matrix(0, nrow = nsim, ncol = maxpq + h)
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
        permanent_component_sim[,i] <- variance_intercept[i] + rho * permanent_component_sim[,i - 1] + phi * (epsilon[,i - 1]^2 - sigma_sqr_sim[,i - 1])
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                transitory_component_sim[,i] <- transitory_component_sim[,i] + alpha[j] * (epsilon[,i - j]^2 - sigma_sqr_sim[,i - j]) + alpha[j] * transitory_component_sim[,i - j]
            }
        }
        if (order[2] > 0) {
            for (j in 1:order[2]) {
                transitory_component_sim[,i] <- transitory_component_sim[,i] + beta[j] * transitory_component_sim[,i - j]
            }
        }
        sigma_sqr_sim[,i] = transitory_component_sim[,i] + permanent_component_sim[,i]
        sigma_sim[,i] <- sqrt(sigma_sqr_sim[,i])
        epsilon[,i] <- z[,i] * sigma_sim[,i]
        series_sim[,i] <- mu + epsilon[,i]
    }
    # return
    sigma <- sigma_sim
    permanent_component <- permanent_component_sim
    transitory_component <- transitory_component_sim
    series <- series_sim
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        permanent_component <- permanent_component[,-seq_len(maxpq + burn), drop = FALSE]
        transitory_component <- transitory_component[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
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
                             vreg = NULL, burn = 0, ...)
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
    extra_args <- list(...)
    h <- h + burn
    parameter <- group <- NULL
    maxpq <- max(object$model$order)
    mu <- object$parmatrix[parameter == "mu"]$value
    omega <- object$parmatrix[parameter == "omega"]$value
    alpha <- object$parmatrix[group == "alpha"]$value
    beta <- object$parmatrix[group == "beta"]$value
    dist <- object$parmatrix[group == "distribution"]$value
    distribution <- object$distribution

    mean_vreg <- 0
    if (!is.null(vreg) & object$vreg$include_vreg) {
        vreg <- as.numeric(vreg)
        if (length(vreg) != h) stop("\nvreg must be a vector of length h + burn.")
        mean_vreg <- mean(vreg)
        variance_intercept <- omega + vreg
        variance_intercept <- c(rep(0, maxpq), variance_intercept)
    } else {
        variance_intercept <- rep(omega, h + maxpq)
    }
    if (object$vreg$multiplicative) variance_intercept <- exp(variance_intercept)

    if (!is.null(innov)) {
        innov <- as.matrix(innov)
        if (NROW(innov) != nsim | NCOL(innov) != h) stop("\ninnov must a matrix of dimensions nsim x (h + burn).")
        z <- cbind(matrix(1, nrow = nsim, ncol = maxpq), innov)
    } else {
        initz <- matrix(1, nrow = nsim, ncol = maxpq)
        z <- matrix(rdist(distribution, h * nsim, mu = 0, sigma = 1, skew = dist[1], shape = dist[2], lambda = dist[3]), nrow = nsim, ncol = h)
        z <- cbind(initz, z)
    }

    if (is.null(var_init)) {
        # set initiale variance close to integrated
        if (object$vreg$multiplicative) {
            numerator <- exp(omega + mean_vreg)
        } else {
            numerator <- omega + mean_vreg
        }
        initv <- numerator/(1 - 0.999)
    } else {
        initv <- var_init
    }

    if (!is.null(innov_init) & maxpq > 0) {
        if (length(innov_init) != maxpq) stop(paste0("\ninnov_init must be of length max(order) : ", maxpq))
        z[,seq_len(maxpq)] <- matrix(innov_init, ncol = maxpq, nrow = nsim, byrow = TRUE)
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

    if (!is.null(extra_args$arch_initial)) {
        init <- extra_args$arch_initial
        if (length(init) != maxpq) init <- rep(init[1], maxpq)
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    } else {
        init <- epsilon[,seq_len(maxpq), drop = FALSE]^2
        init <- matrix(init, ncol = maxpq, nrow = nrow(epsilon), byrow = TRUE)
    }
    simc <- .garchsimvec(z = z, epsilon = epsilon, sigma_sqr_sim = sigma_sqr_sim, variance_intercept = variance_intercept, order = order, init = init, alpha = alpha, beta = beta, mu = mu)
    sigma <- simc$sigma
    series <- simc$series
    if ((maxpq + burn) > 0) {
        sigma <- sigma[,-seq_len(maxpq + burn), drop = FALSE]
        series <- series[,-seq_len(maxpq + burn), drop = FALSE]
    }
    class(sigma) <- "tsmodel.distribution"
    class(series) <- "tsmodel.distribution"
    attr(sigma, "date_class") <- "numeric"
    attr(series, "date_class") <- "numeric"
    out <- list(sigma = sigma, series = series)
    return(out)
}
