test_that("test + operator: multispec",{
    mspec <- global_spec_garch + global_spec_cgarch + global_spec_gjrgarch
    expect_s3_class(mspec, "tsgarch.multispec")
    expect_length(mspec, 3)
    expect_true(attr(mspec, "index_match"))
})

test_that("test as method on estimated objects",{
    L <- list()
    L[[1]] <- global_mod_garch
    L[[2]] <- global_mod_cgarch
    L[[3]] <- global_mod_gjrgarch
    mfit <- to_multi_estimate(L)
    expect_s3_class(mfit, "tsgarch.multi_estimate")
    expect_length(mfit, 3)
    expect_s3_class(mfit[[1]], "tsgarch.estimate")
    expect_s3_class(mfit[[2]], "tsgarch.estimate")
    expect_s3_class(mfit[[3]], "tsgarch.estimate")
})
