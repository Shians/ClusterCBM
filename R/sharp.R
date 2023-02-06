#' Clustering using SHARP
#' @param sce the SingleCellExperiment object
#'
#' @return SingleCellExperiment object
#' @export
cluster_sharp <- function(sce) {
    expr <- logcounts(sce)

    res <- SHARP::SHARP(expr, rN.seed = 1111)
    res <- res$pred_clusters

    colData(sce)$cluster_id <- factor(res)

    return(sce)
}