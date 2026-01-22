load.epic2.manifest <- function() {
    ## download this file https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/InfiniumMethylationEPICv2.0ProductFiles(ZIPFormat).zip
    filename = "MethylationEPIC v2.0 Files/EPIC-8v2-0_A1.csv"

    cat("Reading", basename(filename), "\n")

    require(data.table)
    require(rtracklayer)
    require(GenomicRanges)
    manifest <- fread(filename, skip=7, data.table=F)

    required.columns <- c("SNP_MinorAlleleFrequency","Name","AddressA_ID","AddressB_ID","IlmnID","CHR","MAPINFO")
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

    ## convert coordinates from hg38 to hg19 for consistency with previous manifests
    ## (for some reason illumina only provides hg38)
    chainfile = "hg38ToHg19.over.chain"
    download.file(
      paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/",chainfile,".gz"),
      destfile=paste0(chainfile,".gz"))
    system(paste0("gunzip ", chainfile,".gz"))
    chain = rtracklayer::import.chain("hg38ToHg19.over.chain")
    sites = with(manifest, data.frame(chr=CHR,start=MAPINFO,end=MAPINFO+1))
    rownames(sites) = manifest$IlmnID
    sites = na.omit(sites)
    sites.gr = GenomicRanges::makeGRangesFromDataFrame(sites)  
    sites.hg19 = rtracklayer::liftOver(sites.gr,chain)
    sites.hg19 = as.data.frame(sites.hg19)
    sites.hg19$name = sites.hg19$group_name
    idx = match(sites.hg19$name,manifest$IlmnID)
    manifest$CHR_GRCh38 = manifest$CHR
    manifest$MAPINFO_GRCh38 = manifest$MAPINFO
    manifest$CHR = NA
    manifest$MAPINFO = NA
    manifest$CHR[idx] = as.character(sites.hg19$seqnames)
    manifest$MAPINFO[idx] = sites.hg19$start

    manifest$CHR = sub("^chr","",manifest$CHR)
  
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
