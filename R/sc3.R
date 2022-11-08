# pre-cleaning required for sc3 clustering
#' @importMethodsFrom SingleCellExperiment counts
pre_clean <- function(sce) {
    ave.counts <- rowMeans(counts(sce))
    keep <- ave.counts >= 1
    sce_cleaned <- sce[keep,]
    return(sce_cleaned)
}

#' Clustering using SC3
#' @inheritParams cluster_race_id
#' @param col.sym the column name containing gene symbols
#' @importMethodsFrom SingleCellExperiment rowData colData
#' @importMethodsFrom SummarizedExperiment  rowData<- colData<-
#' @export
cluster_sc3 <- function(sce, col.sym = "Symbol") {
    sce <- pre_clean(sce)
    rowData(sce)$feature_symbol <- rownames(sce)

    sce <- SC3::sc3_prepare(sce, n_cores = 1)
    SingleCellExperiment::isSpike(sce, "ERCC") <- FALSE
    sce <- SC3::sc3_estimate_k(sce)

    k_est <- sce@metadata$sc3$k_estimation
    sce <- SC3::sc3(
        sce,
        ks = k_est,
        biology = FALSE,
        k_estimator = FALSE,
        rand_seed = 0,
        n_cores = 1
    )
    res <- colData(sce)[[glue::glue("sc3_{k_est}_clusters")]]
    res <- factor(res)
    names(res) <- colnames(sce)

    colData(sce)$cluster_id <- res

    return(sce)
}
