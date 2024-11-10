tscombine.garch_modelspec <- function(x, y, ...)
{
    if (!is(y, "tsgarch.spec")) stop("\ny must be of class tsgarch.spec")
    if (is(x, "tsgarch.spec")) {
        anchor_date <- x$target$index
        check_dates <- all.equal(anchor_date, y$target$index)
        index_match <- TRUE
        if (!all(check_dates)) index_match <- FALSE
        out <- list()
        out[[1]] <- x
        out[[2]] <- y
        attr(out,"index_match") <- index_match
        class(out) <- "tsgarch.multispec"
        return(out)
    } else if (is(x,"tsgarch.multispec")) {
        anchor_date <- x[[1]]$target$index
        check_dates <- all.equal(anchor_date, y$target$index)
        index_match <- TRUE
        if (!all(check_dates)) index_match <- FALSE
        x[[length(x) + 1]] <- y
        attr(x,"index_match") <- index_match
        class(x) <- "tsgarch.multispec"
        return(x)
    }
}



#' Combine univariate GARCH specifications into a multi-specification object
#'
#' @param x an object of class \dQuote{tsgarch.spec}
#' @param y an object of class \dQuote{tsgarch.spec}
#' @returns an object of class \dQuote{tsgarch.multispec}
#' @details A simple method for combining multiple specifications into an object
#' which can then be estimated using parallel resources. Note that the returned
#' object is effectively a validated list of specification objects with no names.
#' Names can be assigned post-construction (see example).
#' @method + tsgarch.spec
#' @export
#' @examples
#' library(xts)
#' x <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
#' y <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
#' z <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
#' mspec <- garch_modelspec(x, model = "egarch") +
#' garch_modelspec(y, model = "cgarch") +
#' garch_modelspec(z, model = "aparch")
#' names(mspec) <- c("x", "y", "z")
#' sapply(mspec, function(x) x$model$model)
`+.tsgarch.spec` <- function(x, y)
{
    tscombine.garch_modelspec(x, y)
}

#' @aliases estimate
#' @method estimate tsgarch.multispec
#' @rdname estimate
#' @author Alexios Galanos
#' @export
#'
estimate.tsgarch.multispec <- function(object, solver = "nloptr", control = NULL,
                                       stationarity_constraint = 0.999,
                                       keep_tmb = FALSE, ...)
{
    out <- future_lapply(1:length(object), function(i){
        estimate(object[[i]], solver = solver, control = control, stationarity_constraint = stationarity_constraint,
                 keep_tmb = keep_tmb)
    }, future.packages = "tsgarch", future.seed = TRUE)
    out <- eval(out)
    attr(out, "index_match") <- attr(object, "index_match")
    class(out) <- "tsgarch.multi_estimate"
    return(out)
}



#' Convert a list of tsgarch.estimate objects to a multi_estimate object
#'
#' @param object an list with \dQuote{tsgarch.estimate} objects.
#' @param ... none
#' @returns a validated object of class \dQuote{tsgarch.multi_estimate}.
#' @details
#' This is a convenience method which provides the flexibility
#' to manually estimate each series and then convert the list of these objects
#' to a multi-estimation object, rather than having to use the estimate method
#' on a multi-specification object.
#' @author Alexios Galanos
#' @aliases to_multi_estimate
#' @rdname to_multi_estimate
#' @export
#'
to_multi_estimate <- function(object, ...)
{
    if (!is.list(object)) stop("\nobject must be a list.")
    check_class <- sapply(object, function(x) is(x,'tsgarch.estimate'))
    if (!all(check_class)) {
        invalid_class <- which(!check_class)
        stop(paste0("\nnot all objects in list of class `tsgarch.estimate` :", invalid_class))
    }
    anchor_date <- index(fitted(object[[1]]))
    check_dates <- sapply(2:length(object), function(i) all.equal(anchor_date, index(fitted(object[[i]]))))
    attr(object, "index_match") <- TRUE
    if (!all(check_dates)) attr(object, "index_match") <- FALSE
    class(object) <- "tsgarch.multi_estimate"
    return(object)
}
