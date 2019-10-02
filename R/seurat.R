#' Clustering using Seurat
#' @inheritParams cluster_race_id
#' @export
cluster_seurat <- function(x) {
    counts <- SingleCellExperiment::counts(x)

    sc_data <- Seurat::CreateSeuratObject(counts = counts)
    sc_data <- Seurat::NormalizeData(object = sc_data)
    sc_data@assays$RNA@data <- logcounts(x)
    sc_data <- Seurat::ScaleData(object = sc_data)

    sc_data <- Seurat::FindVariableFeatures(object = sc_data)
    sc_data <- Seurat::RunPCA(
        object = sc_data,
        do.print = FALSE
    )

    sc_data <- Seurat::FindNeighbors(sc_data)
    sc_data <- Seurat::FindClusters(sc_data)

    res <- sc_data@meta.data$seurat_clusters
    res <- factor(res)

    colData(x)$cluster_id <- res

    return(x)
}
