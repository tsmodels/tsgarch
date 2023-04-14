test_that("fcp benchmark",{
    b <- benchmark_fcp()
    expect_true(all(b$lre > 2))
})


# add rugarch benchmark
