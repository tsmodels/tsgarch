.unconditional_garch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$constant_variance * (1 - p)
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}

.unconditional_egarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    unconditional_variance <- exp(numerator/(1 - p))
    return(unconditional_variance)
}

.unconditional_aparch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    delta <- object$parmatrix[parameter == "delta"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    unconditional_variance <- unconditional_variance^(2/delta)
    return(unconditional_variance)
}

.unconditional_gjrgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}

.unconditional_fgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    delta <- object$parmatrix[parameter == "delta"]$value
    p <- persistence(object)
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    unconditional_variance <- unconditional_variance^(2/delta)
    return(unconditional_variance)
}

.unconditional_cgarch <- function(object)
{
    parameter <- group <- NULL
    numerator <- object$parmatrix[parameter == "omega"]$value
    # the permanent component persistence
    p <- persistence(object)[1]
    if (object$spec$model$variance_targeting) {
        numerator <- object$target_omega
    }
    if (object$spec$vreg$include_vreg) {
        xi <- object$parmatrix[group == "xi"]$value
        v <- object$spec$vreg$vreg
        numerator <- numerator + as.numeric(sum(colMeans(v) * xi))
    }
    if (object$spec$vreg$multiplicative) numerator <- exp(numerator)
    unconditional_variance <- numerator/(1 - p)
    return(unconditional_variance)
}
