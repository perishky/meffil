filter.sites <- function(featureset, expr) {
    eval(substitute(meffil:::filter.features(featureset, expr & target == "methylation")))
}

is.autosomal <- function(chr) chr %in% paste0("chr",1:22)
is.chromosome <- function(chr) chr %in% paste0("chr",c(1:22,"X","Y"))

#' @export
meffil.get.typeii.sites <- function(featureset) {
    filter.sites(featureset, type=="ii")
}

#' @export
meffil.get.autosomal.sites <- function(featureset) {
    filter.sites(featureset, is.autosomal(chromosome))
}

#' @export
meffil.get.x.sites <- function(featureset) {
    filter.sites(featureset, chromosome=="chrX")
}

#' @export
meffil.get.y.sites <- function(featureset) {
    filter.sites(featureset, chromosome=="chrY")
}

#' @export
meffil.get.sites <- function(featureset) {
    filter.sites(featureset, !is.na(chromosome))
}

get.quantile.site.subsets <- function(featureset) {
    list(genomic.iG=filter.sites(featureset, type=="i" & meth.dye=="G")$name,
         genomic.iR=filter.sites(featureset, type=="i" & meth.dye=="R")$name,
         genomic.ii=filter.sites(featureset, type=="ii")$name,
         autosomal.iG=filter.sites(featureset, type=="i"&meth.dye=="G"&is.autosomal(chromosome))$name, 
         autosomal.iR=filter.sites(featureset, type=="i"&meth.dye=="R"&is.autosomal(chromosome))$name, 
         autosomal.ii=filter.sites(featureset, type=="ii" & is.autosomal(chromosome))$name, 
         not.y.iG=filter.sites(featureset, type=="i" & meth.dye=="G" & chromosome != "chrY")$name, 
         not.y.iR=filter.sites(featureset, type=="i" & meth.dye=="R" & chromosome != "chrY")$name, 
         not.y.ii=filter.sites(featureset, type=="ii" & chromosome != "chrY")$name, 
         sex=filter.sites(featureset, chromosome %in% c("chrX", "chrY"))$name, 
         chrX=filter.sites(featureset, chromosome=="chrX")$name, 
         chrY=filter.sites(featureset, chromosome=="chrY")$name)
}
         
get.island.site.subsets <- function(featureset) {
    list(island=filter.sites(featureset, relation.to.island=="Island")$name,
         shore=filter.sites(featureset, relation.to.island %in% c("N_Shore","S_Shore"))$name,
         far=filter.sites(featureset, relation.to.island %in% c("N_Shelf","S_Shelf","OpenSea"))$name)
}

is.valid.site.subset <- function(name) {
    name %in% c("genomic.iG", 
                "genomic.iR",
                "genomic.ii",
                "autosomal.iG", 
                "autosomal.iR", 
                "autosomal.ii", 
                "not.y.iG", 
                "not.y.iR", 
                "not.y.ii", 
                "sex", 
                "chry", 
                "chrx") 
}

is.sex.specific.subset <- function(name) {
    name %in% c("genomic.iG",
                "genomic.iR",
                "genomic.ii",
                "not.y.iG",
                "not.y.iR",
                "not.y.ii",
                "sex",
                "chry",
                "chrx")
}

applicable.quantile.site.subsets <- function(sex, both.sexes) {
    if (both.sexes && sex == "M") return(c("autosomal.iG","autosomal.iR","autosomal.ii","sex"))
    if (both.sexes && sex == "F") return(c("autosomal.iG","autosomal.iR","autosomal.ii","chrx","chry"))
    if (!both.sexes && sex == "M") return(c("genomic.iG", "genomic.iR", "genomic.ii"))
    if (!both.sexes && sex == "F") return(c("not.y.iG", "not.y.iR", "not.y.ii","chry"))
    stop("invalid input", "sex =", sex, "both.sexes =", both.sexes)
}
