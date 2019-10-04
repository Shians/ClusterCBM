#' Clustering using Seurat
#' @inheritParams cluster_race_id
#' @export
cluster_seurat <- function(x) {
    counts <- SingleCellExperiment::counts(x)

    pbmc <- Seurat::CreateSeuratObject(
        counts = counts,
        min.cells = 3,
        min.features = 0,
    )

    pbmc <- Seurat::NormalizeData(object = pbmc, display.progress = FALSE)
    pbmc <- Seurat::FindVariableFeatures(
        object = pbmc,
        x.low.cutoff = 0.0125,
        x.high.cutoff = 3,
        y.cutoff = 0.5,
        display.progress = FALSE
    )

    pbmc <- Seurat::ScaleData(
        object = pbmc
    )

    pbmc <- Seurat::RunPCA(
        object = pbmc,
        pcs.compute = num_dim,
        pc.genes = pbmc@var.genes,
        do.print = FALSE
    )

    pbmc <- Seurat::FindNeighbors(pbmc)
    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)

    res <- pbmc@meta.data$seurat_clusters
    res <- factor(res)

    colData(x)$cluster_id <- res

    return(x)
}
