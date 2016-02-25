context("ComputeRP")

test_that("Output of ComputeRP is correct", {
    signif = c("L", "L", "H", "L", "H", "L")
    input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
    rp.list <- ComputeRP(input_file, signif, 100, 1234)
    expect_that(is.data.frame(rp.list), is_true())
    expect_that(ncol(rp.list), equals(4))
    expect_that(is.factor(rp.list[, 1]), is_true())
    expect_that(is.numeric(rp.list[, 2]), is_true())
    expect_that(is.numeric(rp.list[, 3]), is_true())
    expect_that(is.unsorted(rp.list[, 2]), is_false())
    expect_that(is.unsorted(rp.list[, 3]), is_false())
})

test_that("ComputeRP stops if wrong input", {
    input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
    input_dat_RP <- read.table(input_file,header=F, sep='\t',stringsAsFactors = FALSE)
    expect_that(ncol(input_dat_RP)>2, is_true())
    expect_that(is.character(input_dat_RP[, 1]), is_true())
    expect_that(is.numeric(input_dat_RP[, 3]), is_true())
    expect_error(ComputeRP(input_file), "need.*")
    expect_error(ComputeRP(input_file, 2), "signif.*")
    expect_error(ComputeRP(input_file, 2, 100), "signif.*")
    expect_error(ComputeRP(input_file, c("L", "H")), "bad.*")
    expect_error(ComputeRP(input_file, c("L", "H"), 100), "bad.*")
    expect_error(ComputeRP(input_file, c("L", "H", "L", "H", "H", "H"), "equal"), 
        "number.*")
}) 
