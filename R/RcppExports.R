# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_findModes <- function(dataTranspose, HmmVb, nthread) {
    .Call(`_HDclust_rcpp_findModes`, dataTranspose, HmmVb, nthread)
}

rcpp_clust <- function(dataTranspose, HmmVb, rfsClust_, control, nthread) {
    .Call(`_HDclust_rcpp_clust`, dataTranspose, HmmVb, rfsClust_, control, nthread)
}

rcpp_trainHmmVb <- function(dataTranspose, VbStructure, searchControl, trainControl, nthread, VB, HMM, HMMVB, bprint = TRUE) {
    .Call(`_HDclust_rcpp_trainHmmVb`, dataTranspose, VbStructure, searchControl, trainControl, nthread, VB, HMM, HMMVB, bprint)
}

