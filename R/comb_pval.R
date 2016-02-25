#' Convergent evidence based on combined p-values.
#'
#' @description \code{CombP} returns ranks of the genes based on p-values combined using
#'   the famous 'fisher' or 'z-transform' methods.

#' @param file A tab-delimited text file with a minimum of 3 columns. First column should
#'   contain gene names, second column should indicate the evidence type and third column
#'   should contain non-negative numeric values (p-values).
#' @param weight A numeric vector containing weights of evidence types. For example,
#'   sample sizes of various evidence types. If not provided, equal weight is given to all
#'   evidence types.
#' @param method A character string among 'fisher', 'z.transform' or 'logit'.
#' @param na.remove An optional argument, defaults to FALSE. Set this argument to TRUE if
#'   all the genes were not detected across all evidence types.
#' @return If all the inputs are in the correct format as suggested, then the output will
#'   be a dataframe containg genes, their combined p-values and corresponding ranks.
#' @examples
#' cus.weights <- c(100,50,200,300,150,400)
#' input_file_P <- system.file("extdata","CombP_toydata.txt",package="GenRank")
#' CP_ranking <- CombP(input_file_P, method = "fisher", na.remove = TRUE)
#' CP_ranking_z <- CombP(input_file_P, method = "z.transform", na.remove = TRUE, weight = cus.weights)
#' @export
#' @importFrom reshape2 dcast
#' @importFrom survcomp combine.test

# CombP function combines p-values of genes estimated in multiple independent
# studies (evidence types) and returns ranks of genes based on combined p-value.
# When combining p-values of different studies, it is extremely important to ensure
# that all those studies are independently testing the same null hypothesis.

CombP <- function(file, weight, method = c("fisher", "z.transform", "logit"), na.remove = FALSE) {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop("reshape2 package needed for this function to work. Please install it.")
    }
    if (!requireNamespace("survcomp", quietly = TRUE)) {
        stop("survcomp package needed for this function to work. Please install it.")
    }
    requireNamespace("reshape2")
    requireNamespace("survcomp")
    if (missing(file)) {
        stop("No file provided as input")
    }
    if (missing(method)) {
        stop("need to choose a method for combining p-values")
    }
    if (!all(method %in% c("fisher", "z.transform", "logit"))) {
        stop("method can only be one among 'fisher', 'z.transform' and 'logit'")
    }
    file1 <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(file1) < 3) {
        stop("file should contain at least 3 columns")
    }
    if (!all(is.character(file1[, 1]))) {
        stop("the first column of file should only contain gene names")
    }
    if (!all(is.numeric(file1[, 3])) && !all(file1[, 3] > 0)) {
        stop("the third column of file should only contain non-negative numerics")
    }
    # separate genes by evidence type and remove duplicated genes, genes with low p-value
    # are retained in case of duplication
    type.list1 <- split(file1, file1[, 2])
    rm.dup <- function(z){
    z <- z[order(z[,3]),]
    z <- z[!duplicated(z[,1]),]
    }
    sorted.list <- lapply(type.list1, rm.dup)
    sorted.matrix <- do.call("rbind", sorted.list)
    rownames(sorted.matrix) <- NULL
    # generate a dataframe with each row corresponding to a gene's p-values,
    # col-evidence types
    pval.mat <- dcast(sorted.matrix, sorted.matrix[, 1] ~ sorted.matrix[, 2])
    rownames(pval.mat) = pval.mat[, 1]
    pval.mat <- pval.mat[, -1]
    if (any(is.na(pval.mat)) && !na.remove) {
    	stop("NAs were present because genes were not detected across all evidence layers. Set na.remove=TRUE")
    }
    # check whether p-values are available for at least 60% of evidence types
    CheckValid <- function(pvalrow) if (length(which(!is.na(pvalrow)))/length(pvalrow) >
        0.6) {
        TRUE
    } else {
        FALSE
    }
    t.ind <- apply(pval.mat, 1, CheckValid)
    # retain only those genes that were true in the previous check
    pval.filt <- pval.mat[t.ind, ]
    n.del <- nrow(pval.mat) - nrow(pval.filt)
    if (n.del > 0) {
        warning("Genes with no p-values across 60% of evidence types were removed")
    }
    if (missing(weight)) {
        weight <- rep(1, ncol(pval.filt))
    }
    if (nrow(pval.filt) == 0) {
        stop("No genes found with p-values across at least 60% of evidence types")
    }
    # function that can be applied to all rows of pval.filt to combine p-values
    comb.p <- function(matrow) {
        row.weight <- weight[-which(is.na(matrow))]
        combine.test(matrow, method = method, weight = row.weight, na.rm = na.remove)
    }
    comb.pval <- apply(pval.filt, 1, comb.p)
    # return gene list ordered based on rank
    comb.pval <- data.frame(names(comb.pval), comb.pval, rank(comb.pval, ties.method = "min"))
    colnames(comb.pval) = c("Gene", "comb.P", "Rank")
    comb.pval <- comb.pval[order(comb.pval[, 3]), ]
    rownames(comb.pval) <- NULL
    return(comb.pval)
}
