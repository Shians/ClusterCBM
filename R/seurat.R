#' Clustering using Seurat
#' @inheritParams cluster_race_id
#' @importMethodsFrom SingleCellExperiment counts
#' @export
cluster_seurat <- function(x) {
    if (is(x, "SingleCellExperiment")) {
        counts <- SingleCellExperiment::counts(x)
    } else {
        counts <- x
    }

    # heuristically determine number of
    num_dim <- floor(nrow(counts)/5000)
    if (num_dim > 5) {
        num_dim <- 5
    }
    # rough guess of how many PCA components to use
    num_dim <- c(10, 10, 20, 30, 40, 50)[num_dim + 1]
    pbmc <- Seurat::CreateSeuratObject(
        raw.data = counts,
        min.cells = 3,
        min.genes = 0,
        project = "",
        display.progress = FALSE
    )

    pbmc <- Seurat::NormalizeData(object = pbmc, display.progress = FALSE)
    pbmc <- Seurat::FindVariableGenes(
        object = pbmc,
        x.low.cutoff = 0.0125,
        x.high.cutoff = 3,
        y.cutoff = 0.5,
        display.progress = FALSE
    )

    pbmc <- Seurat::ScaleData(
        object = pbmc,
        vars.to.regress = c("nUMI"),
        display.progress = FALSE,
        do.cpp = TRUE
    )

    pbmc <- Seurat::RunPCA(
        object = pbmc,
        pcs.compute = num_dim,
        pc.genes = pbmc@var.genes,
        do.print = FALSE
    )

    pbmc <- Seurat::FindClusters(
        object = pbmc,
        dims.use = 1:num_dim,
        print.output = FALSE
    )

    res <- pbmc@meta.data$res.0.8
    res <- factor(res)
    names(res) <- colnames(counts)

    return(res)
}
