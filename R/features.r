#' @export
meffil.get.features <- function(featureset)
    meffil.featureset(featureset)

filter.features <- function(featureset, expr) {
    features <- meffil.get.features(featureset)
    sat <- eval(substitute(expr), features)
    idx <- which(sat)
    features[idx,]
}



