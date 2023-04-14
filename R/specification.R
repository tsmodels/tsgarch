#' GARCH Model Specification
#'
#' @description Specifies a GARCH model prior to estimation.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines.
#' @param y an xts vector.
#' @param constant whether to estimate a constant (mean) for y,
#' @param model the type of GARCH model. Valid choices are \dQuote{garch} for
#' vanilla GARCH, \dQuote{gjr} for asymmetric GARCH, \dQuote{egarch} for
#' exponential GARCH, \dQuote{aparch} for asymmetric power ARCH,
#' \dQuote{csGARCH} for the component GARCH, \dQuote{igarch} for the integrated
#' GARCH.
#' @param order the (p,q) GARCH order.
#' @param variance_targeting whether to use variance targeting rather than
#' estimating the conditional variance intercept.
#' @param vreg an optional xts matrix of regressors in the conditional variance
#' equation.
#' @param multiplicative whether to exponentiate the contribution of the
#' regressors else will be additive. In the case of the \dQuote{egarch} model,
#' since this is already a multiplicative model, the regressors are additive
#' irrespective of the choice made.
#' @param init the method to use to initialize the recursion of the conditional
#' variance.
#' @param backcast_lambda the decay power for the exponential smoothing used
#' when initializing the recursion using the backcast method.
#' @param sample_n the number of data points to use when initializing the
#' recursion using the sample method.
#' @param distribution a valid distribution from the available
#' re-parameterized distributions of the package.
#' @param ... not used.
#' @return An object of class \dQuote{tsgarch.spec}.
#' @aliases garch_modelspec
#' @rdname garch_modelspec
#' @author Alexios Galanos
#' @export
#'
#'
#'
garch_modelspec <- function(y, model = "garch", constant = FALSE,
                            order = c(1,1), variance_targeting = FALSE,
                            vreg = NULL, multiplicative = FALSE,
                            init = c("unconditional","sample","backcast"),
                            backcast_lambda = 0.7, sample_n = 10,
                            distribution = "norm", ...)

{
    # 1. check and initialize data
    if  (!is.xts(y)) {
        stop("y must be an xts object")
    }
    spec <- initialize_data(y)
    # 2. validate arguments
    model <- match.arg(model[1], choices = valid_garch_models())
    distribution <- match.arg(distribution[1], choices = valid_distributions())
    init <- match.arg(init[1], choices = c("unconditional","sample","backcast"))
    multiplicative <- as.logical(multiplicative)
    constant <- as.logical(constant)
    if (constant) {
        mu <- mean(y, na.rm = TRUE)
    } else {
        mu <- 0.0
    }
    # egarch already in logs
    if (model == "egarch") {
        if (multiplicative) warning("\nmultiplicative not valid for egarch model (already multiplicative due to log specification). Setting to FALSE")
        multiplicative <- FALSE
    }
    # cannot have variance targeting for igarch model
    if (model == "igarch") {
        if (variance_targeting) warning("\nvariance_targeting not possible in igarch model (Inf unconditional variance). Setting to FALSE")
        variance_targeting <- FALSE
    }
    if (model == "ewma") {
        if (variance_targeting) warning("\nvariance_targeting not valid for ewma (omega is zero). Setting to FALSE")
        variance_targeting <- FALSE
    }
    order <- c(order[1], order[2])
    variance_targeting <- as.logical(variance_targeting)

    # 3. check regressors
    if (!is.null(vreg)) {
        vreg <- check_xreg(vreg, index(y))
        if (!multiplicative) {
            if (any(coredata(vreg) < 0)) {
                warning("\nvreg present with negative values and multiplcative = FALSE. Cannot guarantee positivity of variance.")
            }
        }
        if (variance_targeting) {
            warning("\nmultiplicative not available when variance_targeting = TRUE")
            multiplicative <- FALSE
        }
    } else {
        multiplicative <- FALSE
    }
    # 4. populate specification object
    # cmodel: [maxpq, arch, garch, variance_targeting, multiplicative, distribution]
    cmodel <- c(max(order), order[1], order[2], as.integer(variance_targeting),
                as.integer(multiplicative), distribution_class(distribution))
    spec$model$model <- model
    spec$model$order <- order
    spec$model$variance_targeting <- variance_targeting
    spec$model$init <- init
    spec$model$backcast_lambda = backcast_lambda
    spec$model$sample_n <- sample_n
    if (is.null(vreg)) {
        spec$vreg$vreg <- matrix(0, ncol = 1, nrow = NROW(y))
        spec$vreg$include_vreg <- FALSE
        spec$vreg$multiplicative <- multiplicative
    } else {
        spec$vreg$vreg <- coredata(vreg)
        spec$vreg$include_vreg <- TRUE
        spec$vreg$multiplicative <- multiplicative
    }
    spec$distribution <- distribution
    # 5. populate parameters
    parmatrix <- initialize_parameters(model, y, constant = constant,
                                       order = order,
                                       variance_targeting = variance_targeting,
                                       vreg = vreg,
                                       multiplicative = multiplicative,
                                       init = init,
                                       backcast_lambda = backcast_lambda,
                                       sample_n = sample_n,
                                       distribution = distribution)
    if (model == "ewma") spec$model$model <- "igarch"

    spec$parmatrix <- parmatrix
    spec$model_options <- cmodel
    spec$model$constant <- constant
    class(spec) <- "tsgarch.spec"
    return(spec)
}


