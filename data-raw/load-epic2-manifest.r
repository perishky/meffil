load.epic2.manifest <- function() {
    ## https://emea.support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
    ## download this file: https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/MethylationEPIC_v2%20Files.zip

    filename <- "MethylationEPIC v2.0 Files/EPIC-8v2-0_A1.csv"

    cat("Reading", basename(filename), "\n")

    manifest <- read.csv(filename, skip=7, stringsAsFactors=F)

    required.columns <- c("SNP_MinorAlleleFrequency","Name","AddressA_ID","AddressB_ID","IlmnID")
    stopifnot(all(required.columns %in% colnames(manifest)))
    
    ## add snp exclusions
    freq <- manifest$SNP_MinorAlleleFrequency
    freq <- strsplit(freq, ";")
    freq <- lapply(freq, as.numeric)
    L <- sapply(freq,length)
    if (any(L==0))
        freq[which(L==0)] <- 0
    freq <- sapply(freq, max)
    manifest$snp.exclude <- manifest$Name %in% manifest$Name[which(freq > 0.01)]
    ## appears to be SNPs in probe with MAF at least 0.01

    manifest$CHR <- sub("^chr","",manifest$CHR) 

    excluded.addresses <- c(
        "[Controls]",
        "21630339",
        "24669308")
    idx <- which(manifest$AddressA_ID %in% excluded.addresses
                 | manifest$AddressB_ID %in% excluded.addresses
                 | manifest$IlmnID %in% excluded.addresses
                 | manifest$IlmnID %in% excluded.addresses
                 | as.character(manifest$IlmnID) == "[Controls]")
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
