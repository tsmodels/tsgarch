# validation function for currently implemented models
valid_garch_models <- function()
{
    models <- c("garch","gjrgarch","aparch","egarch","fgarch","cgarch","igarch","ewma")
    return(models)
}

# fgarch/avgarch, fgrarch/gjr, fgarch/tgarch, fgarch/ngarch, fgarch/nagarch

# check and process regressors if present for predict
.process_prediction_regressors <- function(old_regressors, new_regressors = NULL, xi = 0, h = 1, include_regressors = FALSE, maxpq = 1, regressor_argument = "newvreg")
{
    if (include_regressors) {
        if (is.null(new_regressors)) stop(paste0("\n",regressor_argument," is NULL but model has a regressor in the variance."))
        if (!is.xts(new_regressors)) new_regressors <- as.matrix(new_regressors)
        if (NROW(new_regressors) != h) stop(paste0("\n",regressor_argument," must have h rows."))
        if (NCOL(new_regressors) != NCOL(old_regressors)) stop(paste0("\n",regressor_argument," must have the same number of columns as regressors in the model."))
        new_v <- rbind(tail(old_regressors,maxpq), coredata(new_regressors))
    } else {
        new_v <- rbind(tail(old_regressors,maxpq), matrix(0, nrow = h, ncol = NCOL(old_regressors)))
    }
    v <- as.numeric(new_v %*% xi)
    return(v)
}


# check and process regressors if present for tsfilter
.process_filter_regressors <- function(old_regressors, new_regressors = NULL, new_index, new_n, include_regressors = FALSE)
{
    n <- length(new_index)
    if (include_regressors) {
        if (is.null(new_regressors)) stop("\nnewvreg is NULL but model has a regressor in the variance.")
        if (NROW(new_regressors) != n) stop("\nnewvreg must have the same number of rows as the new y vector.")
        if (NCOL(new_regressors) != NCOL(old_regressors)) stop("\nnewvreg does not have the same number of columns as vreg in model.")
        # set the index of the new_regressors to that of the new y
        new_regressors <- xts(coredata(new_regressors), new_index)
        v_new <- coredata(.merge_data(old_regressors, new_regressors))
    } else {
        # we do this in case the filtering involves updating old data without appending new data (which is a case
        # not currently supported since it will error out but may be supported in future).
        if (new_n > 0) {
            v_new <- rbind(coredata(old_regressors), matrix(0, nrow = new_n, ncol = NCOL(old_regressors)))
        } else {
            v_new <- coredata(old_regressors)
        }
    }
    return(v_new)
}

# generation of future dates for use in prediction
.forecast_dates <- function(forc_dates = NULL, h = 1, sampling, last_index)
{
    if (is.null(forc_dates)) {
        forc_dates = future_dates(last_index, frequency = sampling, n = h)
    } else {
        if (length(forc_dates) != h) stop("\nforc_dates must be a vector of length h")
        if (any(forc_dates <= last_index)) stop("\nforc_dates must be stricly greater than in-sample dates and times.")
    }
    return(forc_dates)
}

# start variance initialization code options -----------------------------------
initialize_data <- function(y)
{
    n <- NROW(y)
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    sampling <- sampling_frequency(index(y))
    spec <- list()
    spec$target$y_orig <- as.numeric(y)
    spec$target$index <- index(y)
    spec$target$sampling <- sampling
    spec$target$y <- y
    spec$target$good <- good
    return(spec)
}

initialize_variance <- function(y, mu = 0.0, init = c("unconditional","sample","backcast"),
                                backcast_lambda = 0.7, sample_n = 10, delta = 2)
{
    y_demeaned <- y - mu
    v <- switch(init,
                "unconditional" = mean(abs(y_demeaned)^delta),
                "sample" = mean(abs(y_demeaned[1:sample_n])^delta),
                "backcast" = backcast_variance(y_demeaned, backcast_lambda, delta))
    if (delta != 2) {
        v <- v^(2/delta)
    }
    return(v)
}


backcast_variance <- function(y, lambda = 0.7, delta = 2)
{
    # when lambda = 1 == mean(y^2)
    n <- length(y)
    ysqr <- abs(y)^delta
    sigma2 <- mean(ysqr)
    v <- (lambda^n) * sigma2 + (1 - lambda) * sum((lambda^(0:(n - 1)) * ysqr[1:n]))
    return(v)
}

# end variance initialization code options -------------------------------------


# used internally in the estimation code while the TMB object is still alive
# and then can use the estfun method to get the results
score_function <- function(x, env)
{
    # add one call to the fun for models which need to update data
    # (integration done in R for some values since TMB integration remains
    # challenging)
    tmp <- env$fun(x, env)
    - 1.0 * log(env$tmb$report(par = x)$ll_vector)
}


# start imports from the corpcor package ---------------------------------------
.is_positive_definite <- function(x)
{
    x <- as.matrix(x)
    eval <- eigen(x, only.values = TRUE, symmetric = TRUE)$values
    tol <- max(dim(x)) * max(abs(eval)) * .Machine$double.eps
    if (sum(eval > tol) == length(eval)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
.make_positive_definite <- function(x, tol)
{
    x <- as.matrix(x)
    d <- dim(x)[1]
    if (dim(x)[2] != d) stop("Input matrix is not square!")
    es <- eigen(x, symmetric = TRUE)
    esv <- es$values
    if (missing(tol)) tol <- d * max(abs(esv)) * .Machine$double.eps
    delta <- 2 * tol
    tau <- pmax(0, delta - esv)
    dm <- es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(x + dm)
}
# end imports from the corpcor package -----------------------------------------

.lag_vector <- function(x, n_lag = 1, remove_na = FALSE, pad = NA)
{
    # has NAs
    x <- as.matrix(x)
    n <- NROW(x)
    d <- NCOL(x)
    if (d == 1) x <- matrix(x, ncol = 1)
    z <- apply(x, 2, FUN = function(y) .embed_vector(y, n_lag + 1)[,n_lag + 1])
    if (!remove_na) z <- rbind(matrix(pad, ncol = d, nrow = n_lag),z)
    return(z)
}

.embed_vector <- function(x, k, by = 1, ascending = FALSE)
{
    x <- matrix(x, ncol = 1)
    n <- NROW(x)
    s <- seq(1, n - k + 1, by = by)
    lens <- length(s)
    cols <- if (ascending) 1:k else k:1
    return(matrix(x[s + rep(cols, rep(lens,k)) - 1], lens))
}

# for filtering: replace old values with new and append new values
.merge_data <- function(old, new)
{
    m <- NCOL(old)
    index_old <- tail(index(old),1)
    index_new <- index(new)
    if (any(index_new >= index_old)) {
        new <- new[index_new >= index_old]
        merged_data <- merge(old, new)
        inc <- which(is.na(merged_data[,1]))
        coredata(merged_data[inc, c(1:m)]) <- coredata(merged_data[inc, -c(1:m)])
        merged_data <- merged_data[,1:m]
        colnames(merged_data) <- colnames(new)
        return(merged_data)
    } else {
        return(old)
    }
}


.modelspec.tsgarch.estimate <- function(object, ...)
{
    # rebuilt arguments
    args <- list(...)
    if (is.null(args)) {
        spec <- object$spec
        return(spec)
    } else {
        spec <- object$spec
        old_args <- list()
        old_args$y <- spec$target$y
        old_args$model <- spec$model$model
        old_args$order <- spec$model$order
        old_args$constant <- spec$model$constant
        old_args$variance_targeting <- spec$model$variance_targeting
        if (spec$vreg$include_vreg) {
            old_args$vreg <- spec$vreg$include_vreg
        } else {
            old_args$vreg <- NULL
        }
        old_args$multiplicative <- spec$vreg$multiplicative
        old_args$init <- spec$model$init
        old_args$backcast_lambda <- spec$model$backcast_lambda
        old_args$sample_n <- spec$model$sample_n
        old_args$distribution <- spec$distribution
        idx <- chmatch(names(old_args), names(args))
        if (all(is.na(idx))) {
            warnings("\nno matching arguments in ... to garch_modelspec arguments.")
            return(spec)
        } else {
            args_idx <- na.omit(idx)
            old_args_idx <- which(!is.na(idx))
            n <- length(na.omit(idx))
            for (i in 1:n) {
                old_args[[old_args_idx[i]]] <- args[[args_idx[i]]]
            }
            new_spec <- do.call(garch_modelspec, args = old_args, quote = TRUE)
            return(new_spec)
        }
    }
}

.modelspec.tsgarch.spec <- function(object, ...)
{
    # rebuilt arguments
    args <- list(...)
    if (is.null(args)) {
        return(object)
    } else {
        spec <- object
        old_args <- list()
        old_args$y <- spec$target$y
        old_args$model <- spec$model$model
        old_args$order <- spec$model$order
        old_args$constant <- spec$model$constant
        old_args$variance_targeting <- spec$model$variance_targeting
        if (spec$vreg$include_vreg) {
            old_args$vreg <- spec$vreg$include_vreg
        } else {
            old_args$vreg <- NULL
        }
        old_args$multiplicative <- spec$vreg$multiplicative
        old_args$init <- spec$model$init
        old_args$backcast_lambda <- spec$model$backcast_lambda
        old_args$sample_n <- spec$model$sample_n
        old_args$distribution <- spec$distribution
        idx <- chmatch(names(old_args), names(args))
        if (all(is.na(idx))) {
            warnings("\nno matching arguments in ... to garch_modelspec arguments.")
            return(spec)
        } else {
            args_idx <- na.omit(idx)
            old_args_idx <- which(!is.na(idx))
            n <- length(na.omit(idx))
            for (i in 1:n) {
                old_args[[old_args_idx[i]]] <- args[[args_idx[i]]]
            }
            new_spec <- do.call(garch_modelspec, args = old_args, quote = TRUE)
            return(new_spec)
        }
    }
}


# internal function for variance targeting calculation
variance_target <- function(object, ...)
{
    group <- NULL
    pmatrix <- copy(object$parmatrix)
    if (sum(pmatrix$estimate) == 0) {
        # should probably set all of them to 1 for calculations
        pmatrix[group == "alpha", estimate := 1]
        pmatrix[group == "beta", estimate := 1]
    }
    pars <- pmatrix[estimate == 1]$value
    env <- list(distribution = object$distribution, parmatrix = pmatrix, model = object$model$model)
    persist <- .persistence(pars = pars, env)
    sigma_bar <- .sigma_bar(object)
    vreg_bar <- .vreg_bar(object)
    target_omega <- sigma_bar * (1 - persist) - vreg_bar
    return(target_omega)
}


.sigma_bar <- function(object, ...)
{
    group <- parameter <- NULL
    res <- object$target$y_orig - object$parmatrix[parameter == "mu"]$value
    out <- switch(object$model$model,
                  "garch" = mean(res^2),
                  "gjrgarch" = mean(res^2),
                  "aparch" = mean(abs(res)^object$parmatrix[group == "delta"]$value),
                  "fgarch" = mean(abs(res)^object$parmatrix[group == "delta"]$value),
                  "cgarch" = mean(res^2),
                  "egarch" = log(mean(res^2)),
                  "igarch" = mean(res^2))
    return(out)
}

.vreg_bar <- function(object, ...)
{
    group <- NULL
    if (object$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        vreg <- sum(colMeans(object$vreg$vreg) * xi)
    } else {
        vreg <- 0
    }
    return(vreg)
}
