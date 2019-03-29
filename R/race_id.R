#' Clustering using RaceID
#' @param x the SingleCellExperiment or count matrix object to perform clustering on
#'
#' @return vector of cluster identities
#' @export
cluster_race_id <- function(x) {
    if (is(x, "SingleCellExperiment")) {
        counts <- SingleCellExperiment::counts(x)
    } else {
        counts <- x
    }

    sc <- RaceID::SCseq(counts)
    sc <- RaceID::filterdata(sc, mintotal = 1)
    sc <- RaceID::compdist(sc)
    sc <- RaceID::clustexp(sc)
    sc <- RaceID::findoutliers(sc)
    res <- factor(sc@cpart)

    return(res)
}
