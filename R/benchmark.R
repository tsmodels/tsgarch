laurent_benchmark_data <- function()
{
    coefficient <- c(0.04016, 0.04028, 0.15189, 0.46892, 0.84713, 1.33403)
    std_error_h <- c(0.01408, 0.00558, 0.01188, 0.04969, 0.01096, 0.13814)
    out <- data.frame(coefficient, std_error_h)
    rownames(out) <- c("mu","omega","alpha1","gamma1","beta1","delta")
    return(out)
}
#' Laurent APARCH Benchmark
#'
#' @description The APARCH(1,1) benchmark of Laurent (2003).
#' @param control control arguments for the nloptr solver.
#' @returns An object of class \dQuote{benchmark.aparch} which has a \dQuote{as_flextable}
#' method for nice printing of the results.
#' @details The benchmark of Laurent (2003) on the Nikkei daily log returns is
#' based on an APARCH(1,1) model with a constant in the conditional mean
#' equation, and normally distributed errors.
#' @references
#' \insertRef{Laurent2004}{tsgarch}
#' @aliases benchmark_laurent
#' @rdname benchmark_laurent
#' @export
#'
#'
benchmark_laurent <- function(control = nloptr_fast_options())
{

    nikkei <- value <- NULL
    benchmark <- laurent_benchmark_data()
    data("nikkei", envir = environment())
    y <- xts(nikkei[,2], as.Date(nikkei[,1]))
    spec <- garch_modelspec(y, model = "aparch",constant = TRUE, order = c(1,1), init = "unconditional")
    mod <- estimate(spec, control = control)
    coefficient <- unname(coef(mod))
    std_error_h <- sqrt(abs(diag(vcov(mod, type = "H"))))
    out <- data.frame(coefficient, std_error_h)
    rownames(out) <- c("mu","omega","alpha1","gamma1","beta1","delta")
    lre <- do.call(cbind, lapply(1:2, function(i) log_relative_error(x = out[,i], benchmark[,i])))
    rownames(lre) <- rownames(benchmark)
    colnames(lre) <- colnames(benchmark)
    lre <- as.data.frame(lre)
    out <- list(tsgarch = out, benchmark = benchmark, lre = lre,
                symbols = c("\\mu","\\omega","\\alpha_1","\\gamma_1", "\\beta_1","\\delta"),
                headings = c("Coef","Std. Error (H)"))
    class(out) <- "benchmark.laurent"
    return(out)
}

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
#' @returns An object of class \dQuote{benchmark.fcp} which has a \dQuote{as_flextable}
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
    spec <- garch_modelspec(y, model = "garch",constant = TRUE, order = c(1,1), init = "unconditional")
    mod <- estimate(spec, control = control)
    coefficient <- unname(coef(mod))
    std_error_h <- sqrt(abs(diag(vcov(mod, type = "H"))))
    std_error_op <- sqrt(abs(diag(vcov(mod, type = "OP"))))
    std_error_qmle <- sqrt(abs(diag(vcov(mod, type = "QML"))))
    out <- data.frame(coefficient, std_error_h, std_error_op, std_error_qmle)
    rownames(out) <- c("mu","omega","alpha1","beta1")
    lre <- do.call(cbind, lapply(1:4, function(i) log_relative_error(x = out[,i], benchmark[,i])))
    rownames(lre) <- rownames(benchmark)
    colnames(lre) <- colnames(benchmark)
    lre <- as.data.frame(lre)
    out <- list(tsgarch = out, benchmark = benchmark, lre = lre,
                symbols = c("\\mu","\\omega","\\alpha_1","\\beta_1"), headings = c("Coef","Std. Error (H)","Std. Error (OPG)", "Std. Error (QMLE)"))
    class(out) <- "benchmark.fcp"
    return(out)
}



#' Transform an object into flextable
#' @description
#' Transforms a \dQuote{benchmark.fcp} or \dQuote{benchmark.laurent} object into a flextable.
#' @param x an object of class \dQuote{benchmark.fcp} or \dQuote{benchmark.aparch}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... additional arguments passed to flextable method.
#' @returns A flextable object.
#' @importFrom flextable as_flextable
#' @aliases as_flextable.benchmark
#' @method as_flextable benchmark.laurent
#' @rdname as_flextable.benchmark
#' @export
as_flextable.benchmark.laurent <- function(x, digits = 4, ...)
{
    tab <- data.table(" " = c("mu"," ","omega"," ","alpha"," ","gamma"," ","beta"," ","delta"," "),
                      "Estimate.tsgarch" = c(x$tsgarch[1,1],x$lre[1,1],
                                             x$tsgarch[2,1],x$lre[2,1],
                                             x$tsgarch[3,1],x$lre[3,1],
                                             x$tsgarch[4,1],x$lre[4,1],
                                             x$tsgarch[5,1],x$lre[5,1],
                                             x$tsgarch[6,1],x$lre[6,1]
                      ),
                      "Estimate.Laurent" = c(x$benchmark[1,1], as.numeric(NA),
                                         x$benchmark[2,1], as.numeric(NA),
                                         x$benchmark[3,1], as.numeric(NA),
                                         x$benchmark[4,1], as.numeric(NA),
                                         x$benchmark[5,1], as.numeric(NA),
                                         x$benchmark[6,1], as.numeric(NA)
                      ),
                      "Std Error (H).tsgarch" = c(x$tsgarch[1,2],x$lre[1,2],
                                                  x$tsgarch[2,2],x$lre[2,2],
                                                  x$tsgarch[3,2],x$lre[3,2],
                                                  x$tsgarch[4,2],x$lre[4,2],
                                                  x$tsgarch[5,2],x$lre[5,2],
                                                  x$tsgarch[6,2],x$lre[6,2]
                      ),
                      "Std Error (H).Laurent" = c(x$benchmark[1,2], as.numeric(NA),
                                              x$benchmark[2,2], as.numeric(NA),
                                              x$benchmark[3,2], as.numeric(NA),
                                              x$benchmark[4,2], as.numeric(NA),
                                              x$benchmark[5,2], as.numeric(NA),
                                              x$benchmark[6,2], as.numeric(NA)
                      ))
    tab <- flextable(tab) |> separate_header() |> colformat_double(digits = digits) |> italic(i = c(2,4,6,8,10,12)) |> fontsize(i = c(2,4,6,8,10,12), size = 9)
    tab <- tab |> compose(i = 1, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 1, j = 1, as_equation(x$symbols[1]))
    tab <- tab |> compose(i = 3, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 3, j = 1, as_equation(x$symbols[2]))
    tab <- tab |> compose(i = 5, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 5, j = 1, as_equation(x$symbols[3]))
    tab <- tab |> compose(i = 7, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 7, j = 1, as_equation(x$symbols[4]))
    tab <- tab |> compose(i = 9, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 9, j = 1, as_equation(x$symbols[5]))
    tab <- tab |> compose(i = 11, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = 11, j = 1, as_equation(x$symbols[6]))
    tab <- tab |> add_footer_row(colwidths = 5, values = as_paragraph(
        "Notes to table: values in italic are the log relative errors (LRE). The APARCH benchmark is from the paper by Sebastian Laurent, Analytical derivates of the APARCH model, Computational Economics 24 (2004), 51--57."), top = FALSE)
    # replace parameters with symbols
    tab <- tab |> align(i = 1, j = 1:5, align = "justify", part = "footer")
    tab <-  tab |> fontsize(i = 1, size = 9, part = "footer")
    tab <- tab |> set_caption(caption = "APARCH(1,1) Benchmark (Laurent)", html_escape = TRUE)
    return(tab)
}


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
    tab <- tab |> align(i = 1, j = 1:9, align = "justify", part = "footer")
    tab <-  tab |> fontsize(i = 1, size = 9, part = "footer")
    tab <- tab |> set_caption(caption = "GARCH(1,1) Benchmark (FCP)", html_escape = TRUE)
    return(tab)
}


