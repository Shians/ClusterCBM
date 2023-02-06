#' @export
all_cluster_methods <- function() {
    list(
        race_id = cluster_race_id,
        seurat = cluster_seurat,
        tscan = cluster_tscan,
        sharp = cluster_sharp
    )
}
