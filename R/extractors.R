extract_model_values <- function(object, object_type, value_name, ...)
{
    group <- NULL
    parmatrix <- object$parmatrix
    if (object_type == "estimate") {
        x <- object$spec
    } else {
        x <- object
    }
    value <- switch(value_name,
                    "y" = x$target$y,
                    "vreg" = xts(x$vreg$vreg, x$target$index),
                    "mu" = parmatrix[group == "mu"]$value,
                    "omega" = parmatrix[group == "omega"]$value,
                    "phi" = parmatrix[group == "phi"]$value,
                    "rho" = parmatrix[group == "rho"]$value,
                    "alpha" = parmatrix[group == "alpha"]$value,
                    "gamma" = parmatrix[group == "gamma"]$value,
                    "eta" = parmatrix[group == "eta"]$value,
                    "delta" = parmatrix[group == "delta"]$value,
                    "beta" = parmatrix[group == "beta"]$value,
                    "xi" = parmatrix[group == "xi"]$value,
                    "distribution" = parmatrix[group == "distribution"]$value)

    return(value)
}
