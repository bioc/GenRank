# Function that ranks genes within each evidence layer
rankGenes <- function(dataMat) {
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
  dataMat <- dataMat[order(dataMat[, 1]), ]
  rank.Mat <- vector()
  rank.Mat <- cbind(rank.Mat, as.character(dataMat[, 1]))
  rank.Mat <- cbind(rank.Mat, dataMat[, 5])
  return(rank.Mat)
} 
