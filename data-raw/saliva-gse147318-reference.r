#' Defines cell type reference "saliva gse147318"
#' 
#' Adapted for meffil by Haotian Tang
#' Jan 28, 2025

#' Defining a saliva reference panel
#' Based on https://github.com/perishky/meffil/blob/master/data-raw/gse167998-reference.r
#' Ori paper: https://www.tandfonline.com/doi/full/10.1080/15592294.2021.1890874
#'   Middleton, L. Y. M., Dou, J., Fisher, J., Heiss, J. A., Nguyen,
#'   V. K., Just, A. C., … M. Bakulski, K. (2021). Saliva cell type DNA
#'   methylation reference panel for epidemiological studies in
#'   children. Epigenetics, 17(2), 161–177.
#'   https://doi.org/10.1080/15592294.2021.1890874
#' DNA methylation was measured using the Illumina MethylationEPIC BeadChip
#' GSE147318

retrieve.gse147318 <- function(path) {
    source("geo.r")
    dir.create(path)
    ## set max download time
    options(timeout = 1000)
    ## download the data
    filename <- file.path(path, "GSE147318_RAW.tar")
    if (length(list.files(path, pattern = "\\.idat$")) == 0) {
        download.file(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147318/suppl/GSE147318_RAW.tar",
            filename
            )
        cat(date(), "Extracting files from GEO archive.\n")
        system(paste("cd", path, ";", "tar xvf", basename(filename)))
        unlink(filename)
        cat(date(), "Unzipping IDAT files.\n")
        system(paste("cd", path, ";", "gunzip *.idat.gz"))
    } else {
        cat(date(), "IDAT files already exist.\n")
    }
    
    samplesheet <- meffil.create.samplesheet(path)
    samples <- geo.samples("GSE147318")
    names(samples)[names(samples) == "geo_accession"] <- "Sample_Name"
    samples$CellType <- get.characteristic(samples$characteristics_ch1, "cell fraction")
    samples <- merge(samples, samplesheet)
    return(samples)
}

create.saliva.gse147318.reference <- function() {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    if (!requireNamespace("minfi", quietly = TRUE)) {
        BiocManager::install("minfi")
    }
    if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE)) {
        BiocManager::install("IlluminaHumanMethylationEPICmanifest")
    }
    if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
        BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    }

    ## BiocManager::install("minfi")
    require(minfi)
    ## source two scripts from meffil github
    source("minfiEPIC.r")
    `%notin%` <- Negate(`%in%`)
    ## to download saliva data
    ## define directory
    path <- "gse147318"
    #on.exit(unlink(path, recursive = TRUE))
    samples <- retrieve.gse147318(path)
    
    RGset <- read.metharray.exp(targets = samples, force = TRUE)
    
    reference <- preprocessNoob(RGset)
    reference <- mapToGenome(reference)
    
    extracted.data <- extractFromRGSetEPIC(RGset)
    samplesheet <- colData(RGset)
    
    samplesheet <- merge(samplesheet, samples)
    sex <- getSex(reference, cutoff = -3)$predictedSex
    sex <- sign(sex == "F") + 1
    reference <- normalizeFunnormEPIC(
        object = reference, extractedData = extracted.data,
        sex = sex, nPCs = 10, verbose = F
        )
    M <- getMeth(reference)
    U <- getUnmeth(reference)
    
    ## Add cell type reference to meffil
    ## there are 4 different types in samplesheet, and "whole" and "oragene" need to be exclude
    celltypes <- c("CD45pos", "large","whole","oragene")
    cell.types <- celltypes[celltypes %notin% c("whole","oragene")]
    stopifnot(all(cell.types %in% samplesheet$CellType))
    selected <- samplesheet$CellType %in% cell.types
    
    meffil.add.cell.type.reference(
        "saliva gse147318",
        M[, selected], U[, selected],
        cell.types = samplesheet$CellType[selected],
        chip = "epic",
        featureset = "epic",
        description = "Child saliva reference of Middleton et al. Epigenetics 2022",
        verbose = T
        )
}


