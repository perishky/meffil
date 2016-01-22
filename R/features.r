#' Get a list of microarray features from a predefined feature set.
#'
#' @param featureset A name returned by \code{\link{meffil.list.featuresets()}} (Default: "450k").
#' @return A data frame listing all features in the feature set.
#' @export
meffil.get.features <- function(featureset="450k")
    meffil.featureset(featureset)

#' Get a subset of features from a predefined feature set.
#'
#' @param featureset A name returned by \code{\link{meffil.list.featuresets()}}.
#' @param expr An expression using column names in the feature set data frame
#' to select specific features.
#' @return A data frame listing all features satisfying \code{expr}.
#'
#' @examples
#' ## Obtain a list of features measuring methylation using Type II probes.
#' filter.features("450k", target=="methylation" & type == "ii")
#' 
filter.features <- function(featureset, expr) {
    features <- meffil.get.features(featureset)
    sat <- eval(substitute(expr), features)
    idx <- which(sat)
    features[idx,]
}


#' select a compatible featureset for the set of features.
guess.featureset <- function(features) {
    featuresets <- meffil.list.featuresets()
    error.msg <- ""
    sizes <- sapply(featuresets, function(featureset) {
        proposed.features <- meffil.get.features(featureset)$name        
        if (all(features %in% proposed.features))
            length(proposed.features)
        else
            NA
    })
    if (all(is.na(sizes)))
        stop("'features' are not compatible with any available feature set")
    
    featuresets[which.max(sizes)]
}
