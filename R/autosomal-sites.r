#' Restrict to subset of features targetting autosomal CpG sites
#'
#' @param features A vector of feature names. 
#' @return Subset of the input targetting autosomal CpG sites
#' @export
meffil.autosomal.subset <- function(features) {
    all.features <- meffil.all.features()
    features <- all.features[all.features$name %in% features,] 
    features$name[meffil:::is.autosomal(features$chromosome)]
}

