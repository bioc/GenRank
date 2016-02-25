#' Convergent evidence based on rank product method.
#'
#' @description \code{ComputeRP} returns ranks of the genes based on rank product method.
#'
#' @param file A tab-delimited text file with a minimum of 3 columns. First column should
#'   contain gene names, second column should indicate the evidence type and third column
#'   should contain non-negative numeric values (e.g. p-values or effect size).
#' @param signif.type A vector containing letters 'L' or 'H' or both. Length of the vector
#'   should be the same as the number of evidence types. 'L' or 'H' indicate whether the
#'   evidence type contains a low numeric value (e.g. p-value) or a high numeric value
#'   (e.g. effect size).
#' @param n.perm A number indicating number of permutations used to calculate null
#'   density.
#'   Defaults to 100 permutations.
#' @param setseed An optional argument. If provided a numeric value, sampling will be put
#'   in a reproducible state using the setseed value as seed.
#' @return If all the inputs are in the correct format as suggested, then the output will
#'   be a dataframe containg genes, their ranks based on RP and corresponding pfp
#'   (equivalent to FDR).
#' @examples
#' input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
#' signif.val <- c('L','L','H','L','H','L')
#' RP_ranks <- ComputeRP(input_file, signif.type = signif.val)
#' RP_ranks_cust <- ComputeRP(input_file, signif.type = signif.val, n.perm=200, setseed=1234)
#' @export
#' @importFrom matrixStats rowProds

# ComputeRP is a function that ranks genes based on the famous rank-product method
ComputeRP <- function(file, signif.type, n.perm = 100, setseed = NULL) {
    if (missing(file)) {
        stop("No file provided as input")
    }
    if (missing(signif.type)) {
        stop("need to provide signif.type vector")
    }
    if (!all(signif.type %in% c("L", "H"))) {
        stop("signif.type vector can only contain characters 'L' and 'H'")
    }
    if (!all(class(n.perm) == "numeric")) {
        stop("number of permutations (n.perm) should be numeric")
    }
    file1 <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    # Uses a vectorized row product function from matrixStats package
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
        stop("matrixStats package needed for this function to work. Please install it.")
    }
    requireNamespace("matrixStats")
    if (ncol(file1) < 3) {
        stop("file should contain at least 3 columns")
    }
    if (!all(is.character(file1[, 1]))) {
        stop("the first column of file should only contain gene names")
    }
    if (!all(is.numeric(file1[, 3])) && !all(file1[, 3] > 0)) {
        stop("the third column of file should only contain non-negative numerics")
    }
    if (!all(complete.cases(file1))) {
        file1 <- file1[which(complete.cases(file1)), ]
        warning("rows containing NAs were removed")
    }
    n.obs <- table(file1[, 2])
    if (length(signif.type) != length(n.obs)) {
        stop("bad length for signif.type")
    }
    signifType <- rep(signif.type, n.obs)
    file1 <- data.frame(file1, signifType, stringsAsFactors = FALSE)
    colnames(file1)[1:3] = c("Genes", "Origin", "Signif")
    total.genes <- unique(file1[, 1])
    # separately rank genes within each evidence layer add those genes which are missing
    # in each evidence layer to get a rank for each gene
    type.list <- split(file1, file1[, 2])
    ranked.list <- lapply(type.list, function(z) AddGenes(z, total.genes))
    ranks.matrix <- do.call("data.frame", ranked.list)
    rownames(ranks.matrix) <- ranks.matrix[, 1]
    # remove some unnecessary columns that are generated from previous commands
    rm.index <- seq(1, ncol(ranks.matrix) - 1, 2)
    ranks.matrix <- ranks.matrix[, -rm.index]
    # for(i in 1:ncol(ranks.matrix)){ranks.matrix[,i] <-
    # as.numeric(levels(ranks.matrix[,i]))[ranks.matrix[,i]]} convert all columns of
    # this dataframe to numeric
    indx <- sapply(ranks.matrix, is.factor)
    ranks.matrix[indx] <- lapply(ranks.matrix[indx], function(x) as.numeric(levels(x))[x])
    # compute original rank product based on rankings of genes within different evidence
    # layers
    observed.rp <- apply(ranks.matrix, 1, function(rankrow) prod(rankrow)^(1/ncol(ranks.matrix)))
    if (!is.null(setseed)) {
        set.seed(setseed)
    }
    # compare original RP with RP generated based on n.perm permutations compute pfp
    # (equivalent to FDR) and return gene list ordered by top ranks
    top.genes <- ComputePFP(ranks.matrix, n.perm, observed.rp)
    return(top.genes)
}
