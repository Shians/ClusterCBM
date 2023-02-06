cluster_metric <- function(sce, cluster_col, truth_col, method) {
    cluster <- colData(sce)[, cluster_col]
    truth <- colData(sce)[, truth_col]
    igraph::compare(cluster_col, truth_col, method=method)
}

#' Adjusted rand index
#'
#' To use this function in a CellBench workflow, it should be wrapped inside an
#' anonymous function. For example \code{function(x) { cluster_metric_ari(x,
#' "cluster_id", "cell_type") }}
#'
#' @param sce the SingleCellExperiment object
#' @param cluster_col the column name in colData(sce) containing the predicted
#'   clusters
#' @param truth_col the column name in colData(sce) containing the true clusters
#'
#' @return Adjusted Rand index metric for comparing true cluster values and
#'   predicted cluster values.
#' @export
cluster_metric_ari <- function(sce, cluster_col, truth_col) {
    cluster_metric(sce, cluster_col, truth_col, "adjusted.rand")
}

#' Normalized mutual information
#'
#' To use this function in a CellBench workflow, it should be wrapped inside an
#' anonymous function. For example \code{function(x) { cluster_metric_nmi(x,
#' "cluster_id", "cell_type") }}
#'
#' @param sce the SingleCellExperiment object
#' @param cluster_col the column name in colData(sce) containing the predicted
#'   clusters
#' @param truth_col the column name in colData(sce) containing the true clusters
#'
#' @return Normalized mutual information for comparing true cluster values and
#'   predicted cluster values.
#' @export
cluster_metric_nmi <- function(sce, cluster_col, truth_col) {
    cluster_metric(sce, cluster_col, truth_col, "nmi")
}