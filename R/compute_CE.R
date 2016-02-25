#' Convergent Evidence (CE) scores of genes.
#'
#' @description \code{ComputeCE} returns ranks of the genes based on CE scores.
#'
#' @param file A tab-delimited text file with a minimum of 2 columns. First column should
#'   contain gene names and second column should indicate the evidence type.
#' @param PC A character string among 'equal', 'ngene' or 'custom' indicating the prior
#'   credibility.
#' @param cust.weights An optional argument required when the PC='custom'. A numeric
#'   vector containing weights reflecting prior credibility. Should contain as many
#'   weights as the number of evidence types.
#' @return If all the inputs are in the correct format as suggested, then the output will
#'   be a dataframe containing genes and their ranks based on CE scores.
#' @examples
#' input_file <- system.file("extdata","CE_RP_toydata.txt",package="GenRank")
#' CE_ranks <- ComputeCE(input_file,PC = "equal")
#' evid.weight <- c(1,1,0.8,0.8,0.5,1)
#' CE_ranks_cust <- ComputeCE(input_file,PC = "custom", cust.weights = evid.weight)
#' @export

# ComputeCE function computes convergent evidence scores (CE) of genes and returns
# ranks based on CE scores.
ComputeCE <- function(file, PC = c("equal", "ngene", "custom"), cust.weights = NULL) {
    if (missing(file)) {
        stop("No file provided as input")
    }
    if (missing(PC)) {
        stop("need to choose a prior credibility mode (PC)")
    }
    if (!all(PC %in% c("equal", "ngene", "custom"))) {
        stop("PC can only be one among 'equal', 'ngene' and 'custom'")
    }
    PC <- match.arg(PC)
    file1 <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(file1) < 2) {
        stop("file should contain at least 2 columns")
    }
    if (!all(is.character(file1[, 1]))) {
        stop("the first column of file should only contain gene names")
    }
    if (!all(complete.cases(file1))) {
        file1 <- file1[which(complete.cases(file1)), ]
        warning("rows containing NAs were removed")
    }
    # generates a table that shows the presence/absence (binary) of genes across
    # evidence layers
    gene.freq <- table(file1[, 1], file1[, 2])
    # consider genes only once within each evidence layer
    gene.freq[gene.freq > 1] <- 1
    # compute the CE scores depending upon the chosen prior credibility (PC) measure
    switch(PC, equal = {
        CE = apply(gene.freq, 1, mean)
    }, custom = {
        if (is.null(cust.weights)) {
            stop("No custom weights provided")
        }
        n.obs <- table(file1[, 2])
        if (length(cust.weights) != length(n.obs)) {
            stop("bad length for cust.weights")
        }
        if (!all(class(cust.weights) == "numeric")) {
            stop("cust.weights should be numeric")
        }
        weights = cust.weights
        CE = apply(gene.freq, 1, function(geneI) sum(geneI * weights)/sum(weights))
    }, ngene = {
        weights = 1/table(file1[, 2])
        CE = apply(gene.freq, 1, function(geneI) sum(geneI * weights)/sum(weights))
    })
    # sort and rank genes based on observed CE scores
    ce.list <- as.data.frame(CE[order(-CE)])
    rank.list <- rank(-ce.list[, 1], ties.method = "min")
    rank.CE <- cbind(rownames(ce.list), round(ce.list, 3), rank.list)
    rownames(rank.CE) <- NULL
    colnames(rank.CE) <- c("Gene", "CE Score", "Rank")
    return(rank.CE)
}
