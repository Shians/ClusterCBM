#' @export
all_cluster_methods <- function() {
    list(
        race_id = cluster_race_id,
        sc3 = cluster_sc3,
        seurat = cluster_seurat,
        tscan = cluster_tscan
    )
}
