test_that("fcp benchmark",{
    b <- benchmark_fcp()
    expect_true(all(b$lre > 5))
})

test_that("fcp laurent",{
    b <- benchmark_laurent()
    expect_true(all(b$lre$coefficient > 4))
})

