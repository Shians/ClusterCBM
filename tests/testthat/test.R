context("Test all clustering functions")

test_that(
    "Clustering works for SingleCellExperiments", {
    sample_sce_data <- readRDS(cellbench_file("10x_sce_sample.rds"))
    expect_length(cluster_seurat(sample_sce_data), ncol(sample_sce_data))
    expect_length(cluster_sc3(sample_sce_data), ncol(sample_sce_data))
    expect_length(cluster_tscan(sample_sce_data), ncol(sample_sce_data))
    expect_length(cluster_race_id(sample_sce_data), ncol(sample_sce_data))
})

test_that(
    "Clustering works for matrices", {
    library(SingleCellExperiment)
    sample_sce_data <- readRDS(cellbench_file("10x_sce_sample.rds"))
    expect_length(cluster_seurat(counts(sample_sce_data)), ncol(sample_sce_data))
    expect_length(cluster_sc3(counts(sample_sce_data)), ncol(sample_sce_data))
    expect_length(cluster_tscan(counts(sample_sce_data)), ncol(sample_sce_data))
    expect_length(cluster_race_id(counts(sample_sce_data)), ncol(sample_sce_data))
})
