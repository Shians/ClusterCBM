#' Clustering using TSCAN
#' @param sce the SingleCellExperiment object to perform clustering on
#'
#' @return vector of cluster identities
#' @export
cluster_tscan <- function(x) {
    if (is(x, "SingleCellExperiment")) {
        counts <- SingleCellExperiment::counts(x)
    } else {
        counts <- x
    }

    res <- TSCAN::exprmclust(counts)
    res <- res$clusterid
    return(factor(res))
}
