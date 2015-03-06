#library(IlluminaHumanMethylation450kmanifest) ## for getProbeInfo()
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19) ## for probe locations

#' Probe type and location annotation
#'
#' Constructs an data frame annotating probe types
#' for the Infinium HumanMethylation450 BeadChip
#' based on \code{\link[minfi]{getProbeInfo}()}.
#'
#' @param array Microarray identifier (Default: "IlluminaHumanMethylation450k").
#' @param annotation Genomic probe locations annotation (Default: "ilmn12.hg19").
#'
#' @export
meffil.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {
    if (array=="IlluminaHumanMethylation450k" && annotation=="ilmn12.hg19") {
        return(probe.info) ## precomputed, see code in ../data-raw/
    }
    else {
        collect.probe.info(array, annotation)
    }
}

collate.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {
    type1.R <- probe.characteristics("I-Red")
    type1.G <- probe.characteristics("I-Green")
    type2 <- probe.characteristics("II")
    controls <- probe.characteristics("Control")

    msg("reorganizing type information")
    ret <- rbind(data.frame(type="i",target="M", dye="R", address=type1.R$AddressB, name=type1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="M", dye="G", address=type1.G$AddressB, name=type1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii",target="M", dye="G", address=type2$AddressA, name=type2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i",target="U", dye="R", address=type1.R$AddressA, name=type1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="U", dye="G", address=type1.G$AddressA, name=type1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii",target="U", dye="R", address=type2$AddressA, name=type2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressA, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressB, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressA, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressB, name=NA,ext=NA,stringsAsFactors=F),

                 data.frame(type="control",target=controls$Type,dye="R",address=controls$Address, name=NA,ext=controls$ExtendedType,stringsAsFactors=F),
                 data.frame(type="control",target=controls$Type,dye="G",address=controls$Address, name=NA,ext=controls$ExtendedType, stringsAsFactors=F))

    locations <- probe.locations(array, annotation)
    ret <- cbind(ret, locations[match(ret$name, rownames(locations)),])

    ret
}

probe.locations <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {
    annotation <- paste(array, "anno.", annotation, sep="")

    msg("loading probe genomic location annotation", annotation)

    require(annotation,character.only=T)
    data(list=annotation)
    as.data.frame(get(annotation)@data$Locations)
}

probe.characteristics <- function(type) {
    msg("extracting", type)
    minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, type=type)
}
