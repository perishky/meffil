#' Defines cell type reference "blood gse167998"
#' 
#' Salas LA et al. (2022) Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling 

retrieve.gse167998 <- function(path) {
	source("geo.r")

    url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_RAW.tar"
    filename <-  file.path(path, basename(url))
    download.file(url, filename)
    cat(date(), "Extracting files from GEO archive.\n")
    system(paste("cd", path, ";", "tar xvf", basename(filename)))
    unlink(filename)
    cat(date(), "Unzipping IDAT files.\n")
    system(paste("cd", path, ";", "gunzip *.idat.gz"))

	samplesheet <- meffil.create.samplesheet(path)	
    samplesheet <- samplesheet[, c("Sample_Name", "Slide", "Basename")]
    samples <- geo.samples("GSE167998")
    names(samples)[names(samples) == "geo_accession"] <- "Sample_Name"
    names(samples)[names(samples) == "description"] <- "CellType"
    samples <- merge(samples, samplesheet)

    samples
}

create.gse167998.reference <- function() {
	if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!requireNamespace("minfi", quietly = TRUE))
        BiocManager::install("minfi")
    if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
        BiocManager::install("IlluminaHumanMethylationEPICmanifest")
    if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE))
        BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

    require("minfi")
    require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
   	source("minfiEPIC.r")

   	dir.create(path <- "GSE167998")
   	on.exit(unlink(path, recursive=TRUE))
   	samples <- retrieve.gse167998(path)

	RGset <- read.metharray.exp(targets = samples, force=TRUE)

	reference <- preprocessNoob(RGset)
	reference <- mapToGenome(reference)

	extracted.data <- extractFromRGSetEPIC(RGset)
	samplesheet <- colData(RGset)
	samplesheet <- merge(samplesheet, samples)
	sex <- getSex(reference, cutoff = -3)$predictedSex
	sex <- sign(sex=="F")+1
	reference <- normalizeFunnormEPIC(
	    object=reference, extractedData=extracted.data,
	    sex=sex, nPCs=10, verbose=F)
	M <- getMeth(reference)
	U <- getUnmeth(reference)

	## Add cell type reference to meffil
	celltypes <- c("Neu", "Eos", "Bas", "Mono", 
               	   "Bnv", "Bmem", "CD4nv", "CD4mem", 
               	   "Treg", "CD8nv", "CD8mem", "NK", "MIX")
	cell.types <- celltypes[celltypes != "MIX"]
	stopifnot(all(cell.types %in% samplesheet$CellType))
	selected <- samplesheet$CellType %in% cell.types
	meffil.add.cell.type.reference(
	    "blood gse167998",
	    M[,selected], U[,selected],
	    cell.types=samplesheet$CellType[selected],
	    chip="epic",
	    featureset="epic",
	    # specific.sites=IDOLOptimizedCpGs, # use of specific sites would require a licence
	    description="Adult blood reference of Salas et al. Nat Comms 2022",
	    verbose=T)
}

