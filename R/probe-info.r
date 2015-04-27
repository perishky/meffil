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
#' @param verbose If \code{TRUE}, then status messages printed during execution (Default: \code{FALSE}).
#' @return Data frame listing all probes along with annotation information on the micorarray.
#'
#' @export
meffil.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19", verbose=F) {
    msg(verbose=verbose)
    if (array=="IlluminaHumanMethylation450k" && annotation=="ilmn12.hg19") {
        return(probe.info) ## precomputed, see code in ../data-raw/
    }
    else {
        collate.probe.info(array, annotation, verbose=verbose)
    }
}

collate.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19", verbose=F) {

    probe.characteristics <- function(type, verbose=F) {
        msg("extracting", type, verbose=verbose)
        minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, type=type)
    }
    
    type1.R <- probe.characteristics("I-Red", verbose)
    type1.G <- probe.characteristics("I-Green", verbose)
    type2 <- probe.characteristics("II", verbose)
    controls <- probe.characteristics("Control", verbose)
    snps1 <- probe.characteristics("SnpI", verbose)
    snps1.R <- snps1[which(snps1$Color == "Red"),]
    snps1.G <- snps1[which(snps1$Color == "Grn"),]
    snps2 <- probe.characteristics("SnpII", verbose)

    msg("reorganizing type information", verbose=verbose)
    ret <- rbind(data.frame(type="i",target="M", dye="R", address=type1.R$AddressB, name=type1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="M", dye="G", address=type1.G$AddressB, name=type1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii",target="M", dye="G", address=type2$AddressA, name=type2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i-snp",target="MG", dye="R", address=snps1.R$AddressB, name=snps1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i-snp",target="MG", dye="G", address=snps1.G$AddressB, name=snps1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii-snp",target="MG", dye="G", address=snps2$AddressA, name=snps2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i",target="U", dye="R", address=type1.R$AddressA, name=type1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="U", dye="G", address=type1.G$AddressA, name=type1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii",target="U", dye="R", address=type2$AddressA, name=type2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i-snp",target="UG", dye="R", address=snps1.R$AddressA, name=snps1.R$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="i-snp",target="UG", dye="G", address=snps1.G$AddressA, name=snps1.G$Name,ext=NA,stringsAsFactors=F),
                 data.frame(type="ii-snp",target="UG", dye="R", address=snps2$AddressA, name=snps2$Name,ext=NA,stringsAsFactors=F),

                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressA, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressB, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressA, name=NA,ext=NA,stringsAsFactors=F),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressB, name=NA,ext=NA,stringsAsFactors=F),

                 data.frame(type="control",target=controls$Type,dye="R",address=controls$Address, name=NA,ext=controls$ExtendedType,stringsAsFactors=F),
                 data.frame(type="control",target=controls$Type,dye="G",address=controls$Address, name=NA,ext=controls$ExtendedType, stringsAsFactors=F))

    
    annotation <- paste(array, "anno.", annotation, sep="")
    require(annotation,character.only=T)
    data(list=annotation)
    locations <- as.data.frame(get(annotation)@data$Locations)
    islands <- as.data.frame(get(annotation)@data$Islands.UCSC)
    
    ret <- cbind(ret,
                 locations[match(ret$name, rownames(locations)),],
                 islands[match(ret$name, rownames(islands)),])
    ret
}

