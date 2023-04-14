.print_screen <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), table.caption = paste0(toupper(x$model)," Model Summary\n"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters + 1
    coefs <- copy(x$coefficients)
    coef_names <- coefs$term
    coefs <- coefs[,term := NULL]
    coefs <- as.data.frame(coefs)
    rownames(coefs) <- coef_names
    if (!is.null(table.caption)) cat(table.caption)
    cat("\nCoefficients:\n")
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    cat("\nN:", as.integer(x$n_obs))
    cat("\nV(initial):", format(signif(x$initvar, digits = digits)))
    cat(",  ")
    cat("V(unconditional):", format(signif(x$uncvar, digits = digits)))
    cat("\nPersistence:", format(signif(x$persistence, digits = digits)))
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = 2 + digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = 2 + digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = 2 + digits)))
    cat("\n")
    if (x$variance_targeting) cat("\nomega was calculated using variance targeting")
    invisible(x)
}

.print_flextable <- function(x, digits = max(3L, getOption("digits") - 3L),
                             signif.stars = getOption("show.signif.stars"),
                             include.symbols = TRUE, include.equation = TRUE,
                             include.statistics = TRUE,
                             table.caption = paste0(toupper(x$model)," Model Estimation"), ...)
{
    signif <- `Pr(>|t|)` <- NULL
    n <- nrow(x$coefficients)
    cf <- copy(x$coefficients)
    k <- 0
    if (signif.stars) {
        cf[,signif := pvalue_format(`Pr(>|t|)`)]
        out <- flextable(cf) |> set_caption(caption = table.caption) |>
            align(j = "term", align = "left") |>
            align(j = "signif", align = "left") |>
            padding(padding.right = 0, j = "`Pr(>|t|)`", part  = "all") |>
            bold(j = "signif", bold = TRUE) |>
            padding(padding.left = 0, j = "signif", part  = "all") |>
            set_header_labels(term = "", Estimate = "Estimate",
                              `Std. Error` = "Std. Error", `t value` = "t value",
                              `Pr(>|t|)` = "Pr(>|t|)", signif = "" )
        out <- out |> add_footer_lines(values = c("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))
        k <- 1
    } else {
        out <- flextable(cf) |> set_caption(caption = table.caption) |>
            align(j = "term", align = "left") |>
            set_header_labels(term = "", Estimate = "Estimate",
                              `Std. Error` = "Std. Error", `t value` = "t value",
                              `Pr(>|t|)` = "Pr(>|t|)")
    }
    if (include.symbols) {
        for (i in 1:n) {
            out <- compose(out, i = i, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = i,j = 1, as_equation(x$symbol[i]))
        }
    }
    if (include.statistics) {
        out <- out |> add_footer_lines(values = c(paste0("variance targeting: ", x$variance_targeting),
            paste0(sprintf("initialization value: %s",formatC(x$initvar))),
                                                  paste0("LogLik: ", format(x$loglikelihood)),
                                                  sprintf("AIC: %s | BIC: %s", formatC(x$AIC), formatC(x$BIC))))
        k <- k + 4
    }
    if (include.equation) {
        if (k == 0) flag <- TRUE else flag <- FALSE
        out <- out |> add_footer_lines(top = FALSE, values = "Model Equation")
        out <- out |> add_footer_row(top = FALSE, values = " ", colwidths = length(out$col_keys))
        k <- k + 2
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k, j = "term", part = "footer", as_equation(paste0(x$equation[[1]])))
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k + 1, j = "term", part = "footer", as_equation(paste0(x$equation[[2]])))
        k <- k + 2
        out <- out |> add_footer_row(top = FALSE, values = "Persistence (P) and Unconditional Variance Equations", colwidths = length(out$col_keys))
        out <- out |> add_footer_row(top = FALSE, values = " ", colwidths = length(out$col_keys))
        k <- k + 2
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k, j = "term", part = "footer", as_equation(paste0(x$equation[[3]])))
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k + 1, j = "term", part = "footer", as_equation(paste0(x$equation[[4]])))
        k <- k + 2
        if (!flag) out <- out |> hline(part = "footer",i = k - 8, j = 1)
        out <- out |> hline(part = "footer",i = k - 4, j = 1)
        out <- out |> hline(part = "footer",i = k, j = 1)
    }
    out <- colformat_double(out, digits = digits) |> autofit()
    return(out)
}

pvalue_format <- function(x) {
    z <- cut(x, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
    as.character(z)
}
