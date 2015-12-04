load.450k.manifest <- function() {
    filename <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv"

    cat("Reading", filename, "\n")
    system.time(manifest <- read.csv(filename, skip=7, stringsAsFactors=F)) ## 5 minutes

    ## add snp exclusions
    library(openxlsx)
    filename <- "http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx"
    if (!file.exists(basename(filename))) {
        cat("Downloading", filename, "\n")
        download.file(filename, basename(filename))
    }

    cat("Reading sheet 1", "\n")
    cpg.snps <- read.xlsx(basename(filename), sheet=1)

    cat("Reading sheet 2", "\n")
    probe.snps <- read.xlsx(basename(filename), sheet=2)
 
    manifest$snp.exclude <- manifest$Name %in% union(with(cpg.snps, PROBE[which(AF > 0.01)]),
                                                     with(probe.snps, PROBE[which(AF > 0.01)]))


    ## some addresses in the manifest don't appear in idat files
    excluded.addresses <- c("21630339","24669308")
    idx <- which(manifest$AddressA_ID %in% excluded.addresses
                 | manifest$AddressB_ID %in% excluded.addresses
                 | manifest$IlmnID %in% excluded.addresses
                 | manifest$IlmnID %in% excluded.addresses
                 | as.character(manifest$IlmnID) == "[Controls]")
    if (length(idx) > 0)
        manifest <- manifest[-idx,]
        
    manifest
}

