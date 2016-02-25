# Function that add missing genes to each evidence layer
AddGenes <- function(dataMat, allgenes) {
    # ranks genes in each evidence layer ranking (ascending/descending) depends on
    # signifType in col4 of dataMat remove duplicated genes in each evidence layer
    if (dataMat[1, 4] == "L") {
        dataMat <- dataMat[order(dataMat[,3]),]
        dataMat <- dataMat[!duplicated(dataMat[, 1]), ]
        mat.ranks <- rank(dataMat[, 3], ties.method = "min")
    } else if (dataMat[1, 4] == "H") {
        dataMat <- dataMat[order(-dataMat[,3]),]
        dataMat <- dataMat[!duplicated(dataMat[, 1]), ]
        mat.ranks <- rank(-dataMat[, 3], ties.method = "min")
    }
    dataMat <- cbind(dataMat, mat.ranks)
    # missing genes in evidence layer
    Genes <- setdiff(allgenes, dataMat[, 1])
    Origin <- rep(dataMat[1, 2], length(Genes))
    Signif <- rep("NA", length(Genes))
    signifType <- rep("NA", length(Genes))
    # all missing genes given a single modest rank than all possible ranks
    mat.ranks <- rep(length(allgenes) + 1, length(Genes))
    add.genes <- data.frame(Genes, Origin, Signif, signifType, mat.ranks)
    # add missing genes and ranks to existing evidence layer
    dataMat <- rbind(dataMat, add.genes)
    dataMat <- dataMat[order(dataMat[, 1]), ]
    rank.Mat <- vector()
    rank.Mat <- cbind(rank.Mat, as.character(dataMat[, 1]))
    rank.Mat <- cbind(rank.Mat, dataMat[, 5])
    return(rank.Mat)
} 
