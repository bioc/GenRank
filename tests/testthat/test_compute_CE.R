context("ComputeCE")

test_that("Output of ComputeCE is correct", {
    input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
    ce.list <- ComputeCE(input_file, PC = "equal")
    expect_that(is.data.frame(ce.list), is_true())
    expect_that(ncol(ce.list), equals(3))
    expect_that(is.factor(ce.list[, 1]), is_true())
    expect_that(is.numeric(ce.list[, 2]), is_true())
    expect_that(is.integer(ce.list[, 3]), is_true())
    expect_that(is.unsorted(ce.list[, 2][nrow(ce.list):1]), is_false())
    expect_that(is.unsorted(ce.list[, 3]), is_false())
})

test_that("ComputeCE stops if wrong input", {
    input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
    input_dat_RP <- read.table(input_file,header=F, sep='\t',stringsAsFactors = FALSE)
    expect_that(ncol(input_dat_RP)>1, is_true())
    expect_that(is.character(input_dat_RP[, 1]), is_true())
    expect_error(ComputeCE(input_file), "need.*")
    expect_error(ComputeCE(input_file, 2), "PC.*")
    expect_error(ComputeCE("equal"), "need.*")
}) 
