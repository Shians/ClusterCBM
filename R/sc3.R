# pre-cleaning required for sc3 clustering
#' @importMethodsFrom SingleCellExperiment counts
pre_clean <- function(sce) {
    ave.counts <- rowMeans(counts(sce))
    keep <- ave.counts >= 1
    sce1 <- sce[keep,]
    sce1 <- scater::calculateQCMetrics(sce1)
    sce1 <- scater::normalize(sce1)
    return(sce1)
}

#' Clustering using SC3
#' @inheritParams cluster_race_id
#' @param col.sym the column name containing gene symbols
#' @importMethodsFrom SingleCellExperiment rowData colData
#' @importMethodsFrom SummarizedExperiment  rowData<- colData<-
#' @export
cluster_sc3 <- function(sce, col.sym = "Symbol") {
    if (!is(sce, "SingleCellExperiment")) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = sce)
        )
    }

    sce_cleaned <- pre_clean(sce)
    rowData(sce_cleaned)$feature_symbol <- rownames(sce_cleaned)

    sce_cleaned <- SC3::sc3_prepare(sce_cleaned, n_cores = 1)
    SingleCellExperiment::isSpike(sce_cleaned, "ERCC") <- FALSE
    sce_cleaned <- SC3::sc3_estimate_k(sce_cleaned)

    k_est <- sce_cleaned@metadata$sc3$k_estimation
    sce_cleaned <- SC3::sc3(
        sce_cleaned,
        ks = k_est,
        biology = FALSE,
        k_estimator = FALSE,
        rand_seed = 0,
        n_cores = 1
    )
    res <- colData(sce_cleaned)[[glue::glue("sc3_{k_est}_clusters")]]
    res <- factor(res)
    names(res) <- colnames(sce)

    return(res)
}
