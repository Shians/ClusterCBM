#' Clustering using TSCAN
#' @inheritParams cluster_race_id
#' @export
cluster_tscan <- function(sce) {
    expr <- logcounts(sce)

    res <- TSCAN::exprmclust(expr)
    res <- res$clusterid

    colData(sce)$cluster_id <- factor(res)

    return(sce)
}
