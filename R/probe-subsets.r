get.typeii.probes <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$type == "ii")])
}

#' @export
meffil.get.autosomal.sites <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$chr %in% paste0("chr",1:22))])
}

#' @export
meffil.get.x.sites <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$chr == "chrX")])
}

#' @export
meffil.get.y.sites <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$chr == "chrY")])
}

#' @export
meffil.get.sites <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$target %in% c("M","U"))])
}

#' @export
meffil.get.snp.probes <- function() {
    probes <- meffil.probe.info()
    unique(probes$name[which(probes$target %in% c("MG","UG"))])
}

get.quantile.probe.subsets <- function() {
    probes <- meffil.probe.info()
    
    is.iG <- probes$type == "i" & probes$dye == "G"
    is.iR <- probes$type == "i" & probes$dye == "R"
    is.ii <- probes$type == "ii"
    is.genomic <- !is.na(probes$chr)
    is.sex <- is.genomic & probes$chr %in% c("chrX","chrY")
    is.x <- is.genomic & probes$chr == "chrX"
    is.y <- is.genomic & probes$chr == "chrY"
    is.autosomal <- is.genomic & !is.sex
    is.not.y <- is.genomic & probes$chr != "chrY"

    get.probe.subset <- function(in.subset) {
        probes$name[which(probes$target == "M" & in.subset)]
    }

    list(genomic.iG = get.probe.subset(is.iG & is.genomic),
         genomic.iR = get.probe.subset(is.iR & is.genomic),
         genomic.ii = get.probe.subset(is.ii & is.genomic),
         autosomal.iG = get.probe.subset(is.iG & is.autosomal),
         autosomal.iR = get.probe.subset(is.iR & is.autosomal),
         autosomal.ii = get.probe.subset(is.ii & is.autosomal),
         not.y.iG = get.probe.subset(is.iG & is.not.y),
         not.y.iR = get.probe.subset(is.iR & is.not.y),
         not.y.ii = get.probe.subset(is.ii & is.not.y),
         sex = get.probe.subset(is.sex),
         chry = get.probe.subset(is.y),
         chrx = get.probe.subset(is.x))
}

get.island.probe.subsets <- function() {
    probes <- meffil.probe.info()
    
    get.probe.subset <- function(in.subset) {
        probes$name[which(probes$target == "M" & in.subset)]
    }
    ret <- list(island=get.probe.subset(probes$Relation_to_Island == "Island"),
                shore=get.probe.subset(probes$Relation_to_Island %in% c("N_Shore","S_Shore")),
                far=get.probe.subset(probes$Relation_to_Island %in% c("N_Shelf","S_Shelf","OpenSea")))
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

applicable.quantile.probe.subsets <- function(sex, both.sexes) {
    if (both.sexes && sex == "M") return(c("autosomal.iG","autosomal.iR","autosomal.ii","sex"))
    if (both.sexes && sex == "F") return(c("autosomal.iG","autosomal.iR","autosomal.ii","chrx","chry"))
    if (!both.sexes && sex == "M") return(c("genomic.iG", "genomic.iR", "genomic.ii"))
    if (!both.sexes && sex == "F") return(c("not.y.iG", "not.y.iR", "not.y.ii","chry"))
    stop("invalid input", "sex =", sex, "both.sexes =", both.sexes)
}
