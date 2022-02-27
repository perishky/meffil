#' Create a reference for saliva
#' by combining a white blood cell type reference (GEO: GSE35069)
#' and a buccal cell type reference (GEO: GSE48472).
#' 
#' Searched for buccal samples in GEO as follows:
#'  https://gbnci-abcc.ncifcrf.gov/geo/gsm.php
#'  search for GPL acc is equal to GPL13534
#'             Supplementary File contains idat
#'             Characteristics Ch1 contains buccal
retrieve.gse48472 <- function(dir) {
    source("geo.r")
    samples <- geo.samples("GSE48472")
    samples <- samples[grep("buccal", samples$characteristics_ch1),]

    cat("Downloading data ...\n")
    samples$grn.urls <- as.character(samples$supplementary_file)
    samples$red.urls <- as.character(samples$supplementary_file.1)
    for (url in with(samples, c(grn.urls, red.urls))) {
        destfile <- file.path(dir, basename(url))
        download.file(url, destfile=destfile)
        system(paste("gunzip", destfile))
    }

    samples$Sample_Name <- sub("_Grn.idat.gz", "", basename(samples$grn.urls))
    components <- strsplit(samples$Sample_Name, "_")
    samples$Slide <- sapply(components, function(x) x[[2]])
    samples$Array <- sapply(components, function(x) x[[3]])
    samples$Basename <- file.path(dir, samples$Sample_Name)
    samples$CellType <- "Buccal"
    samples$Sex <- "M"
    samples$Sex[which(samples$Sample_Name %in%
                      c("GSM1179538_7310440085_R03C01",
                        "GSM1179540_7310440053_R05C02",
                        "GSM1179539_7310440053_R02C02"))] <- "F"
    samples
}

create.saliva.reference <- function() {
    number.pcs <- 5
    verbose <- T
    
    dir.create(temp.dir <- tempfile(tmpdir="."))
    on.exit(unlink(temp.dir, recursive=TRUE))

    samples.buccal <- retrieve.gse48472(temp.dir)
    samples.blood <- retrieve.gse35069(temp.dir)

    cols <- intersect(colnames(samples.blood), colnames(samples.buccal))

    samplesheet <- rbind(samples.blood[,cols], samples.buccal[,cols])

    ds <- meffil.normalize.dataset(
        samplesheet,
        just.beta=F,
        qc.file="saliva-qc-report.html",
        author="author",
        study="blood and buccal",
        number.pcs=number.pcs,
        norm.file="saliva-normalization-report.html",
        featureset="common",
        verbose=verbose)
    
    samplesheet <- samplesheet[match(colnames(ds$M), samplesheet$Sample_Name),]
    
    cell.types <- c("Buccal","CD4T","CD8T","Mono","Bcell","NK","Gran")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference(
        "saliva gse48472",
        ds$M[,selected], ds$U[,selected],
        cell.types=samplesheet$CellType[selected],
        chip="450k",
        featureset="common",
        description="Saliva reference composed of buccal cell data from Slieker et al. Epigenetics Chromatin 2013 and (blood) immune cell data from Reinius et al. PLoS One 2012",
        verbose=verbose)
}
