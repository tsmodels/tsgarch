fcp_benchmark_data <- function()
{
    coefficient <- c(-0.619041E-2, 0.107613E-1, 0.153134E-0, 0.805974E-0)
    std_error_h <- c(.846212E-2, .285271E-2, .265228E-1, .335527E-1)
    std_error_op <- c(.843359E-2, .132298E-2, .139737E-1, .165604E-1)
    std_error_qmle <- c(.918935E-2, .649319E-2, .535317E-1, .724614E-1)
    out <- data.frame(coefficient, std_error_h, std_error_op, std_error_qmle)
    rownames(out) <- c("mu","omega","alpha1","beta1")
    return(out)
}

log_relative_error <- function(x, b)
{
    -1.0 * log(abs(x - b)/abs(b), base = 10)
}

#' FCP GARCH Benchmark
#'
#' @description The GARCH(1,1) FCP benchmark.
#' @param control control arguments for the nloptr solver.
#' @return An object of class \dQuote{benchmark.fcp} which has a \dQuote{as_flextable}
#' method for nice printing of the results.
#' @details The benchmark of Fiorentini et al. (1996) on the Deutsche Mark British Pound
#' returns is based on a GARCH(1,1) model with a constant in the conditional mean
#' equation, and normally distributed errors.
#' @references
#' \insertRef{Fiorentini1996}{tsgarch}
#' @aliases benchmark_fcp
#' @rdname benchmark_fcp
#' @export
#'
#'
benchmark_fcp <- function(control = nloptr_fast_options())
{
    dmbp <- value <- NULL
    benchmark <- fcp_benchmark_data()
    data("dmbp", envir = environment())
    y <- xts(dmbp[,1], as.Date(1:NROW(dmbp), origin = "1970-01-01"))
    spec <- garch_modelspec(y, model = "garch",constant = TRUE, order = c(1,1))
    mod <- estimate(spec)
    coefficient <- unname(coef(mod))
    std_error_h <- sqrt(abs(diag(vcov(mod, type = "H"))))
    std_error_op <- sqrt(abs(diag(vcov(mod, type = "OP"))))
    std_error_qmle <- sqrt(abs(diag(vcov(mod, type = "QML"))))
    out <- data.frame(coefficient, std_error_h, std_error_op, std_error_qmle)
    rownames(out) <- c("mu","omega","alpha1","beta1")
    benchmark_spec <- garch_modelspec(y, model = "garch",constant = TRUE, order = c(1,1))
    benchmark_spec$parmatrix[estimate == 1, value := benchmark$coefficient]
    benchmark_spec$parmatrix[, estimate := 0]
    # will dispatch to the tsfilter
    filtered_benchmark <- suppressWarnings(estimate(benchmark_spec))
    lre <- do.call(cbind, lapply(1:4, function(i) log_relative_error(x = out[,i], benchmark[,i])))
    rownames(lre) <- rownames(benchmark)
    colnames(lre) <- colnames(benchmark)
    lre <- as.data.frame(lre)
    out <- list(tsgarch = out, benchmark = benchmark, lre = lre, logLik = c(logLik(mod), logLik(filtered_benchmark)),
                symbols = c("\\mu","\\omega","\\alpha_1","\\beta_1"), headings = c("Coef","Std. Error (H)","Std. Error (OPG)", "Std. Error (QMLE)"))
    class(out) <- "benchmark.fcp"
    return(out)
}


#' Transform an object into flextable
#' @description
#' Transforms a \dQuote{benchmark.fcp} object into a flextable.
#' @param x an object of class \dQuote{benchmark.fcp}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... additional arguments passed to flextable method.
#' @return A flextable object.
#' @importFrom flextable as_flextable
#' @aliases as_flextable.benchmark.fcp
#' @method as_flextable benchmark.fcp
#' @rdname as_flextable.benchmark
#' @export
as_flextable.benchmark.fcp <- function(x, digits = 4, ...)
{
    tab <- data.table(" " = c("mu"," ","omega"," ","alpha"," ","beta"," "),
                      "Estimate.tsgarch" = c(x$tsgarch[1,1],x$lre[1,1],
                                             x$tsgarch[2,1],x$lre[2,1],
                                             x$tsgarch[3,1],x$lre[3,1],
                                             x$tsgarch[4,1],x$lre[4,1]
                                             ),
                      "Estimate.FCP" = c(x$benchmark[1,1], as.numeric(NA),
                                         x$benchmark[2,1], as.numeric(NA),
                                         x$benchmark[3,1], as.numeric(NA),
                                         x$benchmark[4,1], as.numeric(NA)
                                         ),
                      "Std Error (H).tsgarch" = c(x$tsgarch[1,2],x$lre[1,2],
                                             x$tsgarch[2,2],x$lre[2,2],
                                             x$tsgarch[3,2],x$lre[3,2],
                                             x$tsgarch[4,2],x$lre[4,2]
                      ),
                      "Std Error (H).FCP" = c(x$benchmark[1,2], as.numeric(NA),
                                         x$benchmark[2,2], as.numeric(NA),
                                         x$benchmark[3,2], as.numeric(NA),
                                         x$benchmark[4,2], as.numeric(NA)
                      ),
                      "Std Error (OPG).tsgarch" = c(x$tsgarch[1,3],x$lre[1,3],
                                                  x$tsgarch[2,3],x$lre[2,3],
                                                  x$tsgarch[3,3],x$lre[3,3],
                                                  x$tsgarch[4,3],x$lre[4,3]
                      ),
                      "Std Error (OPG).FCP" = c(x$benchmark[1,3], as.numeric(NA),
                                              x$benchmark[2,3], as.numeric(NA),
                                              x$benchmark[3,3], as.numeric(NA),
                                              x$benchmark[4,3], as.numeric(NA)
                      ),
                      "Std Error (QML).tsgarch" = c(x$tsgarch[1,4],x$lre[1,4],
                                                    x$tsgarch[2,4],x$lre[2,4],
                                                    x$tsgarch[3,4],x$lre[3,4],
                                                    x$tsgarch[4,4],x$lre[4,4]
                      ),
                      "Std Error (QML).FCP" = c(x$benchmark[1,4], as.numeric(NA),
                                                x$benchmark[2,4], as.numeric(NA),
                                                x$benchmark[3,4], as.numeric(NA),
                                                x$benchmark[4,4], as.numeric(NA)
                      )
                      )
    tab <- flextable(tab) |> separate_header() |> colformat_double(digits = digits) |> italic(i = c(2,4,6,8)) |> fontsize(i = c(2,4,6,8), size = 9)
    tab <- tab |> compose(i = 1, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 1, j = 1, as_equation(x$symbols[1]))
    tab <- tab |> compose(i = 3, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 3, j = 1, as_equation(x$symbols[2]))
    tab <- tab |> compose(i = 5, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 5, j = 1, as_equation(x$symbols[3]))
    tab <- tab |> compose(i = 7, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 7, j = 1, as_equation(x$symbols[4]))
    tab <- tab |> add_footer_row(colwidths = 9, values = as_paragraph(
        "Notes to table: values in italic are the log relative errors (LRE). The FCP benchmark is from the paper by Fiorentini, G., G. Calzolari and L. Panattoni, Analytic Derivatives and the Computation of GARCH Estimates, Journal of Applied Econometrics 11 (1996), 399--417."), top = FALSE)
    # replace parameters with symbols
    tab <-  tab |> fontsize(i = 1, size = 9, part = "footer")
    tab <- tab |> set_caption(caption = "GARCH(1,1) Benchmark (FCP)", html_escape = TRUE)
    return(tab)
}
