#' List of microarrays formats available.
#'
#' By default, there is '450k' and 'epic'.
#' Additions can be made using \code{\link{meffil.add.chip}()}.
#' 
#' @export
meffil.list.chips <- function() {
    ls(probe.globals)
}

#' List of feature sets available.
#'
#' Sets of features for individual platforms (e.g. "450k" for the Illumina HumanMethylation450 Beadchip)
#' as well as for mixed platforms (e.g. "450k:epic:epic2" for
#' combinations of Illumina HumanMethylation450, MethylationEPIC and MethylationEPICv2). 
#'
#' In most cases, a feature corresponds to the two probes
#' from which it's value is derived.  Each CpG represented
#' on the chip for example
#' corresponds to a single feature derived from a probe measuring
#' methylated signal and a second probe measuring unmethylated signal.
#'
#' Each control feature corresponds to a unique control probe.
#' 
#' @export
meffil.list.featuresets <- function() {
    fsnames <- ls(meffil:::featureset.globals)
    fsnames0 <- setdiff(fsnames,"common")
    if (length(fsnames0) > 1) {
        fsnames0 <- unlist(lapply(2:length(fsnames0), function(m) apply(combn(fsnames0,m),2,paste,collapse=":")))
        fsnames <- c(fsnames, fsnames0)
    }
    fsnames
}

#' Obtain a list of features in a feature set.
#'
#' @param featureset Name returned by \code{\link{meffil.list.featuresets}()} (Default: "450k").
#' @return A data frame with one row for each feature.
#' @examples
#' x <- meffil.featureset("450k")
#' 
#' @export
meffil.featureset <- function(featureset="450k") {
   featuresets <- featureset
   is.combo <- !featureset %in% ls(meffil:::featureset.globals)
   if (is.combo)
       featuresets <- unlist(strsplit(featureset, ":"))
   missingsets <- setdiff(featuresets, ls(meffil:::featureset.globals))
   if (length(missingsets) > 0)
       stop(
           paste(missingsets,collapse="/"),
           ifelse(length(missingsets)==1," is ", " are "),
           " not a valid featureset.")
   features <- get(featuresets[1], meffil:::featureset.globals)
   if (is.combo) {
       for (featureset2 in featuresets[-1])
           features <- features[which(features[[featureset2]]),]
   }
   features
}

#' Obtain a list of probes for a given feature set (chip).
#'
#' @param chip Name returned by \code{\link{meffil.list.chips()}}  (Default: "450k").
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{chip}).
#' @return A data frame with one row per probe.  The full set of probes
#' for a chip is returned if \code{chip == featureset}; otherwise,
#' the probes are restricted to those corresponding to features in the feature set.
#' 
#' @export
meffil.probe.info <- function(chip="450k", featureset=chip) {
    probes <- get(chip, probe.globals)
    if (featureset != chip) {
        features <- meffil.featureset(featureset)
        if (!all(features$name %in% probes$name))
            stop(paste("featureset", featureset,
                       "is not compatible with microarray", chip))
        idx <- which(probes$name %in% features$name | probes$target == "OOB")
        probes <- probes[idx,]
    }
    probes
}

#' Add a new chip for analysis.
#'
#' @param name Name of the new chip.
#' @param manifest A data frame obtained by loading the Illumina
#' manifest into R.
#' @return Assuming that \code{manifest} contains a satisfactory
#' set of columns, a new feature set and a new chip is made available.
#' Thus, \code{name} will be added to the vectors returned by
#' \code{\link{meffil.list.featuresets}()} and
#' \code{\link{meffil.list.chips}()}.
#'
#' The manifest must contain the following columns:
#' \itemize{
#' \item{"IlmnID"}{character}
#' \item{"Name"}{character}
#' \item{"AddressA_ID"}{character}
#' \item{"AddressB_ID"}{character}
#' \item{"Infinium_Design_Type"}{values "I","II" or ""}
#' \item{"CHR"}{values "0"-"22", "M", "X" or "Y"}
#' \item{"MAPINFO"}{integer}
#' \item{"AlleleA_ProbeSeq"}{character}
#' \item{UCSC_RefGene_Name}{character}
#' \item{UCSC_RefGene_Accession}{character}
#' \item{UCSC_RefGene_Group}{character}
#' \item{"UCSC_CpG_Islands_Name"}{character}
#' \item{"Relation_to_UCSC_CpG_Island"}{character}
#' \item{"snp.exclude"}{logical}
#' }
#' @export
meffil.add.chip <- function(name, manifest) {
    check.manifest(manifest)
    
    features <- extract.featureset(manifest)
    probes <- extract.probes(manifest)
    meffil.add.featureset(name, features)
    assign(name, probes, probe.globals)

    invisible(manifest)
}

#' Add a feature set.
#'
#' @param name Name of the new feature set.
#' @param features A data frame listing and describing all features.
#' @return Assuming that \code{features} contains a satisfactory
#' set of columns, a new feature set is made available.
#' Thus, \code{name} will be added to the vector returned by
#' \code{\link{meffil.list.featuresets}()}.
#' 
#' The \code{features} data frame must contain the following columns:
#' \itemize{
#' \item{"name"}{character}
#' \item{"target"}{character}
#' \item{"type"}{values "i","ii" or "control"}
#' \item{"chromosome"}{values "chr0"-"chr22", "chrM", "chrX" or "chrY"}
#' \item{"position"}{integer}
#' \item{"gene.symbol"}{character},
#' \item{"gene.accession"}{character},
#' \item{"gene.region"}{character},
#' \item{"cpg.island.name"}{character}
#' \item{"relation.to.island"}{character}
#' \item{"snp.exclude"}{logical}
#' }
#' @export
meffil.add.featureset <- function(name, features) {
    meffil:::check.featureset(features)
    fnames <- with(features, paste(type,target,name))
    for (name2 in meffil.list.featuresets()) {
        features2 <- get(name2, meffil:::featureset.globals)
        fnames2 <- with(features2, paste(type,target,name))
        features2[[name]] <- fnames2 %in% fnames
        features[[name2]] <- fnames %in% fnames2
        assign(name2, features2, meffil:::featureset.globals)
    }
    assign(name, features, meffil:::featureset.globals)
    invisible(features)
}

check.featureset <- function(features) {
    check.data.frame(features,
                     list("type"=c("i","ii","control"),
                          "target"="character",
                          "name"="character",
                          "chromosome"=paste("chr",c(0:22,"X","Y","M"),sep=""),
                          "position"="integer",
                          "gene.symbol"="character",
                          "gene.accession"="character",
                          "gene.region"="character",
                          "cpg.island.name"="character",
                          "relation.to.island"="character",
                          "snp.exclude"="logical"))
}

check.manifest <- function(manifest) {
    stopifnot(is.data.frame(manifest))
    
    check.data.frame(manifest,
                     list("IlmnID"="character",
                          "Name"="character",
                          "AddressA_ID"="character",
                          "AddressB_ID"="character",
                          "Infinium_Design_Type"=c("I","II",""),
                          "CHR"=c(0:22,"X","Y","M",""),
                          "MAPINFO"="integer",
                          "AlleleA_ProbeSeq"="character",
                          "UCSC_RefGene_Name"="character",
                          "UCSC_RefGene_Accession"="character",
                          "UCSC_RefGene_Group"="character",                     
                          "UCSC_CpG_Islands_Name"="character",
                          "Relation_to_UCSC_CpG_Island"="character",
                          "snp.exclude"="logical"))
}

extract.featureset <- function(manifest) {
    manifest$type <- tolower(manifest$Infinium_Design_Type)
    manifest$type[which(!manifest$type %in% c("i","ii"))] <- "control"

    manifest$name <- manifest$Name
    
    manifest$target <- "methylation"
    manifest$target[which(substring(manifest$IlmnID,1,2) == "rs")] <- "snp"

    idx <- which(manifest$type == "control")
    manifest$name[idx] <- manifest$AlleleA_ProbeSeq[idx]
    manifest$target[idx] <- manifest$Name[idx]
    
    ## rename genomic location columns
    manifest$chromosome <- paste("chr", as.character(manifest$CHR), sep="")
    manifest$position <- manifest$MAPINFO
    manifest$chromosome[which(is.na(manifest$position))] <- NA

    ## rename gene columns
    manifest$gene.symbol <- manifest$UCSC_RefGene_Name
    manifest$gene.accession <- manifest$UCSC_RefGene_Accession
    manifest$gene.region <- manifest$UCSC_RefGene_Group
    
    ## rename cpg island columns
    manifest$cpg.island.name <- as.character(manifest$UCSC_CpG_Islands_Name)
    manifest$relation.to.island <- as.character(manifest$Relation_to_UCSC_CpG_Island)
    manifest$relation.to.island[which(manifest$target == "methylation"
                                      & manifest$relation.to.island=="")] <- "OpenSea"

    manifest$meth.dye <- NA
    idx <- which(manifest$type == "i")
    manifest$meth.dye[idx] <- ifelse(manifest$Color_Channel[idx] == "Red","R","G")
    idx <- which(manifest$type == "ii")
    manifest$meth.dye[idx] <- "G"
        
    manifest[,c("type","target","name",
                "chromosome","position", "meth.dye",
                "gene.symbol", "gene.accession","gene.region",
                "cpg.island.name","relation.to.island",
                "snp.exclude")]
}


extract.probes <- function(manifest) {
    type1 <- manifest[which(manifest$Infinium_Design_Type == "I"
                            & substring(manifest$IlmnID,1,2) %in% c("cg","ch")),]
    type1.R <- type1[which(type1$Color_Channel == "Red"),]
    type1.G <- type1[which(type1$Color_Channel == "Grn"),]
    type2 <- manifest[which(manifest$Infinium_Design_Type == "II"
                            & substring(manifest$IlmnID,1,2) %in% c("cg","ch")),]
    controls <- manifest[which(!manifest$Infinium_Design_Type %in% c("I","II")),]
    snps1 <- manifest[which(manifest$Infinium_Design_Type == "I"
                            & substring(manifest$IlmnID, 1,2) == "rs"),]
    snps2 <- manifest[which(manifest$Infinium_Design_Type == "II"
                            & substring(manifest$IlmnID, 1,2) == "rs"),]
    snps1.R <- snps1[which(snps1$Color_Channel == "Red"),]
    snps1.G <- snps1[which(snps1$Color_Channel == "Grn"),]

    probes <- rbind(## methylated channel
                    data.frame(type="i",target="M", dye="R",address=type1.R$AddressB_ID, name=type1.R$Name),
                    data.frame(type="i",target="M", dye="G",address=type1.G$AddressB_ID, name=type1.G$Name),
                    data.frame(type="ii",target="M", dye="G",address=type2$AddressA_ID, name=type2$Name),
                    ## snp 'methylated' channel
                    data.frame(type="i",target="M-snp", dye="R",address=snps1.R$AddressB_ID, name=snps1.R$Name),
                    data.frame(type="i",target="M-snp", dye="G",address=snps1.G$AddressB_ID, name=snps1.G$Name),
                    data.frame(type="ii",target="M-snp", dye="G",address=snps2$AddressA_ID, name=snps2$Name),
                    ## unmethylated channel
                    data.frame(type="i",target="U", dye="R",address=type1.R$AddressA_ID, name=type1.R$Name),
                    data.frame(type="i",target="U", dye="G",address=type1.G$AddressA_ID, name=type1.G$Name),
                    data.frame(type="ii",target="U", dye="R",address=type2$AddressA_ID, name=type2$Name),
                    ## snp 'unmethylated' channel
                    data.frame(type="i",target="U-snp", dye="R",address=snps1.R$AddressA_ID, name=snps1.R$Name),
                    data.frame(type="i",target="U-snp", dye="G",address=snps1.G$AddressA_ID, name=snps1.G$Name),
                    data.frame(type="ii",target="U-snp", dye="R",address=snps2$AddressA_ID, name=snps2$Name),
                    ## out-of-band probes for background adjustment
                    data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressA_ID, name=NA),
                    data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressB_ID, name=NA),
                    data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressA_ID, name=NA),
                    data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressB_ID, name=NA),
                    ## control probes
                    data.frame(type="control", target=controls$Name, dye="R",address=controls$IlmnID, name=controls$AlleleA_ProbeSeq),
                    data.frame(type="control",target=controls$Name,dye="G",address=controls$IlmnID, name=controls$AlleleA_ProbeSeq))

    for (i in 1:ncol(probes))
        probes[,i] <- as.character(probes[,i])
    
    probes
}


is.compatible.chip <- function(featureset, chip) {
    ret <- FALSE
    try({
        probes <- meffil.probe.info(chip, featureset)
        ret <- TRUE
    }, silent=TRUE)
    ret
}

# select a compatible chip the \code{object}
# which may be an 'rg' object (see \code{\link{read.rg}()}
# or a matrix (typically beta or methylation or unmethylation matrices)
# with row names corresponding to feature/probe names.
guess.chip <- function(object, chip=NA) {
    stopifnot(is.character(chip) || is.na(chip))
    
    chips <- chip
    if (is.na(chip))
        chips <- meffil.list.chips()
    error.msg <- ""
    if (is.rg(object)) {
        for (chip in chips) {
            probes <- meffil.probe.info(chip)
            if (all(probes$address %in% c(rownames(object$G), rownames(object$R))))
                return(chip)
        }
        error.msg <- paste("IDAT files (", object$basename, ")", sep="")
    } else if (is.matrix(object)) {
        for (chip in chips) {
            features <- meffil.featureset(chip)
            if (all(rownames(object) %in% features$name))
                return(chip)
        }
        error.msg <- "Matrix"
    }
    error.msg <- paste(error.msg, "is not compatible with the following chip(s):",
                       paste(chips, collapse=", "))
    stop(error.msg)
}
