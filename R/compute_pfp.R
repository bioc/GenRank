# Function that generates permutation based RP, computes PFP and returns top genes
# performs permutation n.perm times computes pfp (refer to original rankprod
# package/article/wiki)
ComputePFP <- function(rankmatrix, n.perm = 100, orig.rp) {
    # pval <- vector() length(pval) = nrow(rankmatrix)
    pval <- numeric(nrow(rankmatrix))
    all.PermRP <- replicate(n.perm, GenPerm(rankmatrix))
    for (i in 1:nrow(all.PermRP)) {
        pval[i] = length(all.PermRP[i, ][all.PermRP[i, ] <= orig.rp[i]])/n.perm
    }
    # return gene list ordered based on rank
    top.list <- data.frame(names(orig.rp), round(orig.rp, 3), rank(orig.rp, ties.method = "min"), 
        round(pval/rank(orig.rp), 3))
    colnames(top.list) = c("Gene", "RP", "Ranks", "pfp")
    top.list <- top.list[order(top.list[, 3]), ]
    rownames(top.list) <- NULL
    return(top.list)
}
 
