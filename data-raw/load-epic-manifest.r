load.epic.manifest <- function() {
    ## filename <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v1-0-b1-manifest-file-csv.zip"
    ## filename <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v1-0-b2-manifest-file-csv.zip"
    ## filename <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v1-0-b3-manifest-file-csv.zip"
    filename <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip"

    on.exit(unlink(basename(filename)))
    cat("Downloading", filename, "\n")
    download.file(filename, basename(filename))

    cat("Reading", basename(filename), "\n")

    manifest.filename <- unzip(basename(filename),list=T)$Name[1]  
    system.time(manifest <- read.csv(unz(basename(filename), manifest.filename),
                                     skip=7, stringsAsFactors=F)) ## 5 minutes
    
    ## add snp exclusions
    manifest$snp.exclude <- manifest$Name %in% with(manifest, Name[which(as.numeric(as.character(SNP_MinorAlleleFrequency)) > 0.01)])

    ## remove leading 0's from probes addresses
    manifest$AddressA_ID <- sub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- sub("^0*", "", manifest$AddressB_ID)

    ## some addresses in the manifest don't appear in idat files
    ## excluded.addresses <- readLines("epic-excludes.txt")
    ## idx <- which(manifest$AddressA_ID %in% excluded.addresses
    ##              | manifest$AddressB_ID %in% excluded.addresses
    ##              | manifest$IlmnID %in% excluded.addresses
    ##              | manifest$IlmnID %in% excluded.addresses
    ##              | as.character(manifest$IlmnID) == "[Controls]")
    idx <- which(as.character(manifest$IlmnID) == "[Controls]")
    if (length(idx) > 0)
         manifest <- manifest[-idx,]
    
    manifest
}
