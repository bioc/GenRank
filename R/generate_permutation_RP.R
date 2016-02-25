# Given a rank matrix, this function permutes the gene labels once Subsequentlt this
# function computes the RP of permuted rank matrix Returns the resultant RP vector
# (length=number of genes)
GenPerm <- function(rankmatrix) {
    perm.rankmatrix <- apply(rankmatrix, 2, function(x) sample(x, length(x), FALSE))
    # perm.rp <- apply(perm.rankmatrix,1,function(rankrow)
    # prod(rankrow)^(1/ncol(perm.rankmatrix)))
    perm.rp <- (rowProds(perm.rankmatrix))^(1/ncol(perm.rankmatrix))
    return(perm.rp)
}
 
