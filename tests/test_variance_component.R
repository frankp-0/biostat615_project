library(data.table)
library(testthat)
source("../code/variance_component.R")

                                        # Test 1
test_that("Test 1", {
    df <- fread("vc_1.txt")
    n <- nrow(Sg)
    ind <- seq(1, 2 * n, 2)
    Sg <- as.matrix(fread("../data/sg.txt.gz"))
    Ce <- as.matrix(fread("../data/re.txt.gz"))
    sqrt_Sg_ginv <- sqrt_ginv(Sg)
    res <- run_vc_optimizer(df, sqrt_Sg_ginv, Ce, ind)
    expect_equal(res$pleio_stat, 407067.231365, tolerance = 1e-9)
    })
