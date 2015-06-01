# meffil
Efficient algorithms for analyzing DNA methylation data.

The `minfi` version consists of a single function call that requires
access to all data. The `meffil` version splits the method up into
several functions in order to allow parallelization and to reduce the
amount of data that needs to be loaded.


## Overview

	# Load meffil and set how many cores to use for parallelization
	library(meffil)
	options(mc.cores=6)

	# Generate samplesheet
	samplesheet <- meffil.create.samplesheet(path_to_idat_files)

	# Or read in samplesheet
	samplesheet <- meffil.read.samplesheet(path_to_idat_files)

	# Background and dye bias correction, sexprediction, cell counts estimates
	qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE)

    # Obtain genotypes for comparison with those measured on the microarray
	genotypes <- meffil.extract.genotypes(plink.files)

	# Generate QC report
	qc.summary <- meffil.qc.summary(qc.objects, genotypes=genotypes)
	meffil.qc.report(qc.summary, output.file="qc-report.html")

	# Remove outlier samples if necessary
	qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

	# Perform quantile normalization
	norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

	# Generate normalized probe values
	norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)

	# Generate normalization report
	norm.summary <- meffil.normalization.summary(norm.beta, norm.objects)
	meffil.normalization.report(norm.summary, output.file="normalization-report.html")

## More info

Get some data:

	dir.create(path <- "~/data/test_meffil", recursive=TRUE)

	if (length(list.files(path, "*.idat$")) == 0) {
	  filename <-  file.path(path, "gse55491.tar")
	  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
	  cat(date(), "Extracting files from GEO archive.\n")
	  system(paste("cd", path, ";", "tar xvf", basename(filename)))
	  unlink(filename)
	  cat(date(), "Unzipping IDAT files.\n")
	  system(paste("cd", path, ";", "gunzip *.idat.gz"))
	}

Load up `meffil`

	library(devtools)
	install_github("perishky/meffil")
	library(meffil)
	options(mc.cores=16)

These data don't actually have a samplesheet, so we can generate one from the `.idat` basenames:

	samplesheet <- meffil.create.samplesheet(path)

The function creates the following necessary columns:

- Sample_Name
- Sex (possible values "M"/"F"/NA)
- Basename

And it also tries to parse the basenames to guess if the Sentrix plate and positions are present. At this point it is worthwhile to manually modify the `samplesheet` data.frame to replace the actual sample IDs in the `Sample_Name` column if necessary, and to add the sex values to the `Sex` column. Don't change these column names though.

Perform the background correction, dye bias correction, sex prediction and cell count estimates:

	qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE)

A list of available cell type references can be obtained as follows:

	meffil.get.cell.type.references()

New references can be created from a dataset using meffil.create.cell.type.reference().

Obtain the matrix of genotypes for comparison with those measured on the microarray.
If such a matrix is available (rows = SNPs, columns = samples), then the following steps
can be omitted.  Otherwise, it is possible to obtain the matrix from a PLINK
dataset as follows:

    writeLines(meffil.get.snp.probes(), con="snp-names.txt")
    command shell > plink -bfile dataset --extract snp-names.txt --recodeA --out genotypes.raw --noweb
    filenames <- "genotypes.raw"
    genotypes <- meffil.extract.genotypes(filenames)

We can now summarise the QC analysis of the raw data

	qc.summary <- meffil.qc.summary(qc.objects, genotypes=genotypes)

and generate a report:
	
	meffil.qc.report(qc.summary, output.file="qc-report.html")

This creates the file "qc-report.html" in the current work directory.
Should open up in your web browser.

You can remove bad samples prior to performing quantile normalization:

	qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

And now perform quantile normalization:

	norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

Finally, the `beta` values can be generated, whilst removing CpGs that were found to be dodgy in the QC analysis:

	norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.report$bad.cpgs$name)

A summary report of the normalization performance can also be generated:

	norm.summary <- meffil.normalization.summary(norm.beta, norm.objects)
	meffil.normalization.report(norm.summary, output.file="normalization-report.html")




