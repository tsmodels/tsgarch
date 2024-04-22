#' Walk Forward Rolling Backtest
#'
#' @description Generates an expanding window walk forward backtest with option
#' for rolling the forecast by filtering (see details).
#' @param object an object of class \dQuote{tsgarch.spec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param estimate_every number of periods at which the model is re-estimated
#' (defaults to 1).
#' @param rolling this indicates whether forecasts are made only on the estimation
#' date (FALSE) or whether to filter the data 1 period at a time and forecast
#' from the filtered data (TRUE).
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... not currently used.
#' @returns A list which includes a data.table having the following columns:
#' \itemize{
#' \item estimation_date: the date at which the model was estimated.
#' \item convergence: whether both kkt1 and kkt2 were TRUE (
#' Kuhn Karush Tucker conditions) from the kktchk function in optimx.
#' \item filter_date: the date on which a prediction was generated. For rolling
#' prediction this means that an estimated model was filtered for new data prior
#' to re-predicting.
#' \item horizon: the forecast horizon of the prediction.
#' \item size: the length of the data used in estimation.
#' \item forecast_date: the date corresponding to the forecast.
#' \item mu: the conditional mean prediction.
#' \item sigma: the conditional volatility prediction.
#' \item skew: the distribution skew parameter (non-time varying hence constant
#' across each estimation window).
#' \item shape: the distribution shape parameter (non-time varying hence constant
#' across each estimation window).
#' \item shape: the distribution lambda parameter (non-time varying hence constant
#' across each estimation window).
#' \item actual: the actual observation corresponding to the forecast date.
#' }
#' Additional slots in the list include the distribution used and other information
#' relating to the backtest setup.
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @details The rolling option allows to re-estimate the data every n periods
#' whilst filtering the data 1-step ahead between re-estimation dates so that overlapping
#' forecasts are generated.
#' @aliases tsbacktest
#' @method tsbacktest tsgarch.spec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.tsgarch.spec <- function(object, start = floor(length(object$target$y_orig))/2,
                                    end = length(object$target$y_orig), h = 1,
                                    estimate_every = 1, rolling = FALSE, trace = FALSE, ...)
{

    # return mu, sigma, skew, shape, lambda
    parameter <- b <- forecast_dates <- NULL
    data <- xts(object$target$y_orig, object$target$index)
    if (object$vreg$include_vreg) {
        use_vreg <- TRUE
        vreg <- xts(object$vreg$vreg, object$target$index)
    } else {
        use_vreg <- FALSE
        vreg <- NULL
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    idates <- index(data)
    # check for ending period to avoid overlap
    index_table <- rbindlist(lapply(1:length(seqdates), function(i) {
        if (i == length(seqdates)) {
            s <- index(data[paste0(seqdates[i],"/",max(idates[end]))])
            s <- s[-length(s)]
        } else {
            s <- index(data[paste0(seqdates[i],"/",seqdates[i + 1])])
            s <- s[-length(s)]
        }
        if (length(s) == 0) return(NULL)
        rbindlist(lapply(1:length(s), function(j){
            idx <- which(idates == s[j])
            maxh <- min(idx + h, NROW(data[1:end]))
            tmp <- na.omit(index(data[idx:maxh]))
            data.table(estimate_date = seqdates[i], filter_date = tmp[1], forecast_dates = tmp[-1], h = length(tmp[-1]))
        }))
    }))
    max_index <- max(idates[end])
    index_table <- index_table[forecast_dates <= max_index]
    index_table <- split(index_table, by = "estimate_date")
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    i <- 1
    b %<-% future_lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        if (use_vreg) {
            vreg_train <- vreg[index(y_train)]
        } else {
            vreg_train <- NULL
        }
        spec <- garch_modelspec(y_train, constant = object$model$constant, order = object$model$order,
                                variance_targeting = object$model$variance_targeting, vreg = vreg_train,
                                multiplicative = object$vreg$multiplicative, init = object$model$init,
                                backcast_lambda = object$model$backcast_lambda, sample_n = object$model$sample_n,
                                distribution = object$distribution)
        mod <- try(estimate(spec), silent = TRUE)
        model_coef <- coef(mod)
        skew <- mod$parmatrix[parameter == "skew"]$value
        shape <- mod$parmatrix[parameter == "shape"]$value
        lambda <- mod$parmatrix[parameter == "lambda"]$value
        converge <- mod$conditions$kkt1 * mod$conditions$kkt2
        L <- index_table[[i]]
        M <- split(L, by = "filter_date")
        rolls <- length(M)
        # first index is the estimation so we start with the prediction
        # vreg
        if (use_vreg) {
            vreg_test <- vreg[M[[1]]$forecast_dates]
        } else {
            vreg_test <- NULL
        }
        # actual
        P <- vector(mode = "list", length = rolls)
        tmp <- predict(mod, h = M[[1]]$h[1], newvreg = vreg_test, forc_dates = M[[1]]$forecast_date, nsim = 0)
        if (rolls > 1 & rolling) {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            P[[1]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                       "convergence" = rep(converge, nh),
                       "filter_date" = M[[1]]$filter_date,
                       "horizon" =  1:M[[1]]$h[1],
                       "size" = rep(nrow(y_train), nh),
                       "forecast_date" =  M[[1]]$forecast_dates,
                       "mu" = as.numeric(tmp$mean),
                       "sigma" = as.numeric(tmp$sigma),
                       "skew" = rep(skew, nh),
                       "shape" = rep(shape, nh),
                       "lambda" = rep(lambda, nh),
                       "actual" = as.numeric(y_test))
            for (j in 2:rolls) {
                roll_dates <- M[[j]]$filter_date[1]
                updated_y <- data[paste0(seqdates[i],"/",roll_dates)][-1,]
                if (use_vreg) {
                    newvreg <- vreg[paste0(seqdates[i],"/",roll_dates)][-1,]
                } else {
                    newvreg <- NULL
                }
                f <- tsfilter(mod, y = updated_y, newvreg = newvreg)
                if (use_vreg) {
                    vreg_test <- vreg[M[[j]]$forecast_dates]
                } else {
                    vreg_test <- NULL
                }
                tmp <- predict(f, h = M[[j]]$h[1], newvreg = vreg_test, forc_dates = M[[j]]$forecast_dates, nsim = 0)
                y_test <- data[M[[j]]$forecast_dates]
                nh <- length(M[[j]]$h)
                P[[j]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                  "convergence" = rep(converge, nh),
                                  "filter_date" = M[[j]]$filter_date,
                                  "horizon" =  1:M[[j]]$h[1],
                                  "size" = rep(nrow(y_train), nh),
                                  "forecast_date" =  M[[j]]$forecast_dates,
                                  "mu" = as.numeric(tmp$mean),
                                  "sigma" = as.numeric(tmp$sigma),
                                  "skew" = rep(skew, nh),
                                  "shape" = rep(shape, nh),
                                  "lambda" = rep(lambda, nh),
                                  "actual" = as.numeric(y_test))
            }
            out <- rbindlist(P)
        } else {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            out <- data.table("estimation_date" = rep(seqdates[i], nh),
                              "convergence" = rep(converge, nh),
                              "filter_date" = M[[1]]$filter_date,
                              "horizon" =  1:M[[1]]$h[1],
                              "size" = rep(nrow(y_train), nh),
                              "forecast_date" =  M[[1]]$forecast_dates,
                              "mu" = as.numeric(tmp$mean),
                              "sigma" = as.numeric(tmp$sigma),
                              "skew" = rep(skew, nh),
                              "shape" = rep(shape, nh),
                              "lambda" = rep(lambda, nh),
                              "actual" = as.numeric(y_test))
        }
        return(out)
    }, future.packages = c("tsmethods","tsgarch","xts","data.table"), future.stdout	= FALSE, future.seed = FALSE)
    b <- eval(b)
    b <- rbindlist(b)
    out <- list(table = b, distribution = object$distribution, h = h, estimate_every = estimate_every, rolling = rolling)
    return(out)
}

