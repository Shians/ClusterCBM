#' Clustering using RaceID
#' @param sce the SingleCellExperiment object
#'
#' @return SingleCellExperiment object
#' @export
cluster_race_id <- function(sce) {
    sc <- RaceID::SCseq(logcounts(sce))
    sc <- RaceID::filterdata(sc, mintotal = 1)
    sc <- RaceID::compdist(sc)
    sc <- RaceID::clustexp(sc)
    sc <- RaceID::findoutliers(sc)

    res <- factor(sc@cpart)

    colData(sce)$cluster_id <- res

    return(sce)
}
