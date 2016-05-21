filter.sites <- function(featureset, expr) {
    eval(substitute(meffil:::filter.features(featureset, expr & target == "methylation")))$name
}

is.autosomal <- function(chr) chr %in% paste0("chr",1:22)
is.chromosome <- function(chr) chr %in% paste0("chr",c(1:22,"X","Y"))

#' Get names of CpG sites corresponding to Infinium Type II probes
#' in the feature set.
#' 
#' @export
meffil.get.typeii.sites <- function(featureset="450k") {
    filter.sites(featureset, type=="ii")
}

#' Get names of autosomal CpG sites in the feature set.
#' 
#' @export
meffil.get.autosomal.sites <- function(featureset="450k") {
    filter.sites(featureset, is.autosomal(chromosome))
}


#' Get names of chromosome X CpG sites in the feature set.
#' @export
meffil.get.x.sites <- function(featureset="450k") {
    filter.sites(featureset, chromosome=="chrX")
}

#' Get names of chromosome Y CpG sites in the feature set.
#' @export
meffil.get.y.sites <- function(featureset="450k") {
    filter.sites(featureset, chromosome=="chrY")
}

#' Get names of all CpG sites in the feature set.
#' @export
meffil.get.sites <- function(featureset="450k") {
    filter.sites(featureset, !is.na(chromosome))
}

get.quantile.site.subsets <- function(featureset) {
    list(genomic.iG=filter.sites(featureset, type=="i" & meth.dye=="G"),
         genomic.iR=filter.sites(featureset, type=="i" & meth.dye=="R"),
         genomic.ii=filter.sites(featureset, type=="ii"),
         autosomal.iG=filter.sites(featureset, type=="i"&meth.dye=="G"&is.autosomal(chromosome)), 
         autosomal.iR=filter.sites(featureset, type=="i"&meth.dye=="R"&is.autosomal(chromosome)), 
         autosomal.ii=filter.sites(featureset, type=="ii" & is.autosomal(chromosome)), 
         not.y.iG=filter.sites(featureset, type=="i" & meth.dye=="G" & chromosome != "chrY"), 
         not.y.iR=filter.sites(featureset, type=="i" & meth.dye=="R" & chromosome != "chrY"), 
         not.y.ii=filter.sites(featureset, type=="ii" & chromosome != "chrY"), 
         sex=filter.sites(featureset, chromosome %in% c("chrX", "chrY")), 
         chrx=filter.sites(featureset, chromosome=="chrX"), 
         chry=filter.sites(featureset, chromosome=="chrY"))
}
         
get.island.site.subsets <- function(featureset) {
    list(island=filter.sites(featureset, relation.to.island=="Island"),
         shore=filter.sites(featureset, relation.to.island %in% c("N_Shore","S_Shore")),
         far=filter.sites(featureset, relation.to.island %in% c("N_Shelf","S_Shelf","OpenSea")))
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
