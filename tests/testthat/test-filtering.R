
test_that("filter from estimate (ghst) append: check results",{
    filt <- tsfilter(global_mod_ghst, y = y[1801,1], newvreg = y[1801,2])
    expect_equal(filt$sigma[1:1800], global_mod_ghst$sigma)
    expect_equal(residuals(global_mod_ghst, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1801)
    expect_length(filt$spec$target$y_orig, 1801)
    expect_length(as.numeric(filt$spec$target$y), 1801)
    expect_equal(filt$nobs, 1801)
})

test_that("filter from spec (ghst) no append: check results",{
    global_spec_ghst$parmatrix <- copy(global_mod_ghst$parmatrix)
    filt <- tsfilter(global_spec_ghst)
    expect_equal(filt$sigma[1:1800], global_mod_ghst$sigma)
    expect_equal(residuals(global_mod_ghst, standardize = TRUE),
                 residuals(filt, standardize = T)[1:1800], ignore_attr = TRUE)
    expect_length(filt$sigma, 1800)
    expect_length(filt$spec$target$y_orig, 1800)
    expect_length(as.numeric(filt$spec$target$y), 1800)
    expect_equal(filt$nobs, 1800)
    expect_equal(filt$loglik, global_mod_ghst$loglik)
})
