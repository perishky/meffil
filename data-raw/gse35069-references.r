#' Defines cell type references "blood gse35069", "blood gse35069 complete"
#' and "blood gse35069 chen".
#'
#' for estimating blood cell counts.
#' Both use methylation profiles from Reinius et al. 2012 (PMID 25424692)
#' for purified blood cell types.
#' The first is based on
#' six cell types: CD4T, CD8T, Mono, Bcell, NK, Gran.
#' The second is based on 
#' the same cell types but with Gran replaced by Neu and Eos.
#'
#' An epigenome-wide association study of total serum IgE in Hispanic children.
#' Chen W, Wang T, Pino-Yanes M, Forno E, Liang L, Yan Q, Hu D, Weeks DE,
#' Baccarelli A, Acosta-Perez E, Eng C, Han YY, Boutaoui N, Laprise C,
#' Davies GA, Hopkin JM, Moffatt MF, Cookson WO, Canino G, Burchard EG, Celedón JC.
#' J Allergy Clin Immunol. 2017 Jan 6. pii:
#' S0091-6749(16)32546-5. PMID: 28069425


retrieve.gse35069 <- function(dir) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(dir)

    #' GEO does not have the IDAT files
    #' (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069);
    #' however, the `FlowSorted.Blood.450k` R package
    #' contains the following code to download these
    #' files from another website as well as the sample
    #' information (see `getKereData.R` script).
    
    cat("Downloading data ...\n")
    download.file("http://publications.scilifelab.se/f742bd8a08ad42c3921ccaaaf0e3997a/file/sample_sheet_IDAT.csv", destfile="samplesheet.csv", quiet=T)
    download.file("http://publications.scilifelab.se/f742bd8a08ad42c3921ccaaaf0e3997a/file/idat_files.zip", destfile="idat-files.zip",quiet=T)

    cat("Unzipping data ...\n")
    unzip("idat-files.zip")

    samples <- read.csv("samplesheet.csv", stringsAsFactors = TRUE)
    names(samples)[names(samples) == "Sample"] <- "Sample_Name"
    names(samples)[names(samples) == "Chip.ID"] <- "Slide"
    samples$Array <- paste0(samples$Chip.Row.Pos, samples$Chip.CO.position)
    samples$Chip.CO.position <- samples$Chip.Row.Pos <- NULL
    samples$Array <- gsub("O", "0", samples$Array)
    samples$Basename <- file.path(dir, paste0(samples$Slide, "_", samples$Array))
    samples$CellType <- samples$Type
    samples$CellType <- sub("^CD14\\+ Monocytes", "Mono", samples$CellType)
    samples$CellType <- sub("^CD19\\+ B-cells", "Bcell", samples$CellType)
    samples$CellType <- sub("^CD4\\+ T-cells", "CD4T", samples$CellType)
    samples$CellType <- sub("^CD56\\+ NK-cells", "NK", samples$CellType)
    samples$CellType <- sub("^CD8\\+ T-cells", "CD8T", samples$CellType)
    samples$CellType <- sub("^Granulocytes", "Gran", samples$CellType)
    samples$CellType <- sub("^Whole blood", "WBC", samples$CellType)
    samples$CellType <- sub("^Neutrophils", "Neu", samples$CellType)
    samples$CellType <- sub("^Eosinophils", "Eos", samples$CellType)
    samples$Sex <- "M"

    ## check.samplesheet(samples)
    
    samples
}

create.gse35069.references <- function() {
    number.pcs <- 5
    verbose <- T
    
    dir.create(temp.dir <- tempfile(tmpdir="."))
    on.exit(unlink(temp.dir, recursive=TRUE))
    
    samplesheet <- retrieve.gse35069(temp.dir)

    ds <- meffil.normalize.dataset(samplesheet,
                                   just.beta=F,
                                   qc.file="gse35069-qc-report.html",
                                   author="Reinius, et al.",
                                   study="Purified blood cell type methylation (GEO:GSE35069)",
                                   number.pcs=number.pcs,
                                   norm.file="gse35069-normalization-report.html",
                                   featureset="common",
                                   verbose=verbose)

    samplesheet <- samplesheet[match(colnames(ds$M), samplesheet$Sample_Name),]
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Neu","Eos")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference("blood gse35069 complete",
                                   ds$M[,selected], ds$U[,selected],
                                   cell.types=samplesheet$CellType[selected],
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Gran")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference("blood gse35069",
                                   ds$M[,selected], ds$U[,selected],
                                   cell.types=samplesheet$CellType[selected],
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)


    cpg.sites <- read.csv("chen-table-e2.csv",stringsAsFactors=F)    

    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Neu","Eos")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference("blood gse35069 chen",
                                   ds$M[,selected], ds$U[,selected],
                                   cell.types=samplesheet$CellType[selected],
                                   chip="450k",
                                   featureset="common",
                                   specific.sites=cpg.sites$cpg,
                                   verbose=verbose)    
}




