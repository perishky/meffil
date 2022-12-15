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
    freq <- sapply(freq, max)
    manifest$snp.exclude <- manifest$Name %in% manifest$Name[which(freq > 0.01)]
    ## appears to be SNPs in probe with MAF at least 0.01

    manifest$CHR <- sub("^chr","",manifest$CHR) 
    
    manifest
}
