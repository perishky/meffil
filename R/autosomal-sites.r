autosomal.sites <- function(beta) {
    if (is.matrix(beta))
        all.sites <- rownames(beta)
    else {
        all.sites <- retrieve.gds.cpg.sites(beta)   
    }    
    featureset <- meffil:::guess.featureset(all.sites)
    autosomal.sites <- meffil.get.autosomal.sites(featureset)
    intersect(autosomal.sites, all.sites)
}


