load.mouse.manifest <- function() {
    ## https://emea.support.illumina.com/array/array_kits/infinium-mouse-methylation-beadchip-kit/downloads.html
    ## https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/infinium-mouse-methylation-manifest-file-csv.zip
    ## https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv
    
    filename <- "MouseMethylation-12v1-0_A2.csv"
    annotation.filename <- "MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv"
    cat("Reading", basename(filename), "\n")

    manifest <- read.csv(filename, skip=7, stringsAsFactors=F)
    annotation <- read.csv(annotation.filename,stringsAsFactors=F)
    
    required.columns <- c("Name","AddressA_ID","AddressB_ID","IlmnID")
    stopifnot(all(required.columns %in% colnames(manifest)))
    
    manifest$snp.exclude <- F
    
    idx <- match(manifest$Name, annotation$name)
    manifest$UCSC_RefGene_Name <- annotation$Gene[idx]
    manifest$UCSC_RefGene_Accession <- ""
    manifest$UCSC_RefGene_Group <- annotation$Feature[idx]
    manifest$UCSC_CpG_Islands_Name <- manifest$CpG_Island
    manifest$Relation_to_UCSC_CpG_Island <- ""

    manifest$Infinium_Design_Type <- ifelse(manifest$Infinium_Design_Type==1,"I","II")
    
    manifest$CHR <- sub("^chr","",manifest$CHR)
    manifest$CHR[which(manifest$CHR=="MT")] <- "M"

    idx <- which(as.character(manifest$IlmnID) == "[Controls]")
    if (length(idx) > 0)
         manifest <- manifest[-idx,]

    ## a few thousand cytosines are measured multiple times on the microarray
    ## using the same probe sequences each time
    ## we'll make all available but append the duplicated probes with the Illumina ID
    ## which includes a duplication number
    probe.names <- manifest$Name[grep("^(cg|ch|rs)", manifest$Name)]
    probe.freq <- table(probe.names)
    duplicate.probes <- names(which(probe.freq > 1))
    dup.idx <- which(manifest$Name %in% duplicate.probes)
    dup.idx <- setdiff(dup.idx, match(duplicate.probes, manifest$Name)) ## omit the first time a probe appears
    manifest$Name[dup.idx] <- manifest$IlmnID[dup.idx]
    ## check
    stopifnot(all(duplicate.probes %in% manifest$Name)) ## all duplicated probes appear >= 1 time
    stopifnot(all(table(manifest$Name[grep("^(cg|ch|rs)", manifest$Name)]))==1) ## all probe names appear 1 time
    
    manifest
}
