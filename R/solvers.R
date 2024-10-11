#' Default options for nloptr solver
#'
#' @param trace equivalent to option \dQuote{print_level} for nloptr. High values results in
#' more details.
#' @param xtol_rel Stop when an optimization step (or an estimate of the optimum) changes
#' every parameter by less than xtol_rel multiplied by the absolute value of the parameter.
#' @param xtol_abs is a vector of length n (the number of elements in x) giving the
#' tolerances: stop when an optimization step (or an estimate of the optimum)
#' changes every parameter x(i) by less than xtol_abs(i).
#' @param maxeval top when the number of function evaluations exceeds maxeval.
#' This is not a strict maximum: the number of function evaluations may exceed
#' maxeval slightly, depending upon the algorithm
#' @return A list with options which can be passed to the solver.
#' @details These as just a set of pre-created defaults which work well, particularly
#' the \dQuote{nloptr_fast_options} which uses an SQP solver. nloptr has many other
#' solvers and combination of solvers which can be used. However, keep in mind that the solver
#' must accept analytic derivatives as well as nonlinear inequality constraints.
#' @export nloptr_fast_options
#' @aliases nloptr_options
#' @rdname nloptr_options
#'
nloptr_fast_options <- function(trace = FALSE, xtol_rel = 1e-14, maxeval = 1000, xtol_abs = 1e-12)
{
    opt <- list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_SLSQP", xtol_rel = xtol_rel, maxeval = maxeval, xtol_abs = xtol_abs, check_derivatives = FALSE)
    return(opt)
}

#' @export nloptr_global_options
#' @aliases nloptr_options
#' @rdname nloptr_options
#'
nloptr_global_options <- function(trace = FALSE, xtol_rel = 1e-14, maxeval = 1000, xtol_abs = 1e-12)
{
    opt <- list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_AUGLAG", xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval, check_derivatives = FALSE,
                local_opts = list(algorithm = "NLOPT_LD_MMA", maxeval = 500, xtol_rel = 1e-12))
    return(opt)
}

# function adapted from the optimx::kktchk
solver_conditions <- function(pars, fn, gr, hess, lower, upper, env)
{
    kkttol <- 0.001
    kkt2tol <- 1e-06
    kkt1 <- NA
    kkt2 <- NA
    npar <- length(pars)
    nbm <- 0
    fval <- fn(pars, env)
    ngr <- gr(pars, env)
    nHes <- hess(pars, env)
    pHes <- nHes
    gmax <- max(abs(ngr))
    kkt1 <- (gmax <= kkttol * (1 + abs(fval)))
    phev <- try(eigen(pHes)$values, silent = TRUE)
    if (!inherits(phev, "try-error")) {
        negeig <- (phev[npar] <= (-1) * kkt2tol * (1 + abs(fval)))
        evratio <- phev[npar - nbm]/phev[1]
        kkt2 <- (evratio > kkt2tol) && (!negeig)
        ans <- list(evratio, kkt1, kkt2)
        names(ans) <- c("evratio", "kkt1", "kkt2")
        return(ans)
    }
    else {
        evratio <- NA
        ans <- list(evratio, kkt1, kkt2)
        names(ans) <- c("evratio", "kkt1", "kkt2")
        return(ans)
    }
}



