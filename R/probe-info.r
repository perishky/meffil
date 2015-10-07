#library(IlluminaHumanMethylation450kmanifest) ## for getProbeInfo()
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19) ## for probe locations

#' Probe type and location annotation
#'
#' @return The data frame being used to annotate
#' the Infinium HumanMethylation450 BeadChip.
#' The default uses \code{\link[minfi]{getProbeInfo}()}
#' and \code{\link[IlluminaHumanMethylation450k]{IlluminaHumanMethylation450kmanifest}}
#' extensively.  Probe locations are obtained using
#' \code{\link[IlluminaHumanMethylation450kanno.ilmn12.hg19]{IlluminaHumanMethylation450kanno.ilmn12.hg19}}.
#' 
#' @export
meffil.probe.info <- function() {
    if (!exists("probe.info", pkg.globals))
        stop("probe.info object not created, see data-raw/globals.r")
    get("probe.info", pkg.globals)
}

#' Change probe type and location annotation
#'
#' This function will fail if expected columns are missing.
#' See the existing annotation by calling \code{\link{meffil.probe.info}()}.
#' 
#' @export
meffil.set.probe.info <- function(probe.info) {
    stopifnot(is.data.frame(probe.info))

    if (exists("probe.info", pkg.globals)) {
        columns <- colnames(get("probe.info", pkg.globals))
        if (length(setdiff(columns, colnames(probe.info))) > 0)
            stop(paste("Annotation is missing several columns:", paste(columns, collapse=",")))
    }
    
    assign("probe.info", probe.info, pkg.globals)
}



collate.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19", verbose=F) {

    probe.characteristics <- function(type, verbose=F) {
        # msg("extracting", type, verbose=verbose)
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
    snpinfo <- as.data.frame(get(annotation)@data$SNPs.137CommonSingle)
    snpinfo$snp_exclude <- with(snpinfo, (CpG_maf > 0.01 | Probe_maf > 0.01))
    snpinfo$snp_exclude[is.na(snpinfo$snp_exclude)] <- FALSE
    snpinfo <- subset(snpinfo, select=c(snp_exclude))
    
    ret <- cbind(ret,
                 locations[match(ret$name, rownames(locations)),],
                 islands[match(ret$name, rownames(islands)),],
                 snpinfo[match(ret$name, rownames(snpinfo)),]
            )
    ret
}

