#' @export
meffil.list.architectures <- function() {
    ls(probe.globals)
}

#' @export
meffil.list.featuresets <- function() {
    ls(featureset.globals)
}

#' @export
meffil.featureset <- function(featureset) {
   stopifnot(featureset %in% meffil.list.featuresets())
   get(featureset, featureset.globals)
}

#' @export
meffil.probe.info <- function(featureset, architecture) {
    if (missing(architecture))
        architecture <- featureset
    
    probes <- get(architecture, probe.globals)
    if (featureset != architecture) {
        features <- meffil.featureset(featureset)
        if (!all(features$name %in% probes$name))
            stop(paste("featureset", featureset,
                       "is not compatible with microarray", architecture))
        idx <- which(probes$name %in% features$name | probes$target == "OOB")
        probes <- probes[idx,]
    }
    probes
}

#' @export
meffil.add.manifest <- function(name, manifest) {
    check.manifest(manifest)
    
    features <- extract.featureset(manifest)
    probes <- extract.probes(manifest)

    assign(name, features, featureset.globals)
    assign(name, probes, probe.globals)

    return(TRUE)
}

#' @export
meffil.add.featureset <- function(name, features) {
    check.featureset(features)
    assign(name, features, featureset.globals)
}

check.featureset <- function(features) {
    check.data.frame(features,
                     list("type"=c("i","ii","control"),
                          "target"="character",
                          "name"="character",
                          "chromosome"=paste("chr",c(1:22,"X","Y"),sep=""),
                          "position"="integer",
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
                          "CHR"=c(1:22,"X","Y",""),
                          "MAPINFO"="integer",
                          "AlleleA_ProbeSeq"="character",
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


is.compatible.architecture <- function(featureset, architecture) {
    ret <- FALSE
    try({
        probes <- meffil.probe.info(featureset, architecture)
        ret <- TRUE
    }, silent=TRUE)
    ret
}

guess.architecture <- function(object) {
    if (is.rg(object)) {
        for (architecture in meffil.list.architectures()) {
            probes <- meffil.probe.info(architecture, architecture)
            if (all(probes$address %in% c(rownames(object$G), rownames(object$R))))
                return(architecture)
        }
    } else if (is.matrix(object)) {
        for (architecture in meffil.list.architectures()) {
            features <- meffil.featureset(architecture)
            if (all(rownames(object) %in% features$name))
                return(architecture)
        }
    }
    return(FALSE)
}
