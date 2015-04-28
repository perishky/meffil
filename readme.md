# meffil
Efficient algorithms for analyzing DNA methylation data.

The `minfi` version consists of a single function call that requires access to all data. The `meffil` version splits the method up into several functions in order to allow parallelization and to reduce the amount of data that needs to be loaded.


## Overview

	# Load meffil and set how many cores to use for parallelization
	library(meffil)
	options(mc.cores=6)

	# Generate samplesheet
	samplesheet <- meffil.create.samplesheet(path_to_idat_files)

	# Or read in samplesheet
	samplesheet <- meffil.read.samplesheet(path_to_idat_files)

	# Background and dye bias correction, sexprediction
	qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

	# Generate QC report
	qc.summary <- meffil.qc.summary(qc.objects)
	meffil.qc.report(qc.summary)

	# Remove outlier samples if necessary
	qc.objects <- meffil.remove.ids(qc.objects, qc.summary$bad.ids)

	# Perform quantile normalization
	qc.quantiles <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

	# Generate normalized probe values
	normalized.beta <- meffil.normalize.samples(qc.quantiles, cpglist.remove=qc.report$bad.cpgs$name)

	# Generate normalization report
	normalization.summary <- meffil.normalization.summary(normalized.beta, qc.quantiles)
	meffil.normalization.report(normalization.summary)


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
	github_install("perishky/meffil@restructure1")
	library(meffil)
	options(mc.cores=16)

These data don't actually have a samplesheet, so we can generate one from the `.idat` basenames:

	samplesheet <- meffil.create.samplesheet(path)

The function creates the following necessary columns:

- Sample_Name
- Sex
- Basenames

And it also tries to parse the basenames to guess if the Sentrix plate and positions are present. At this point it is worthwhile to manually modify the `samplesheet` data.frame to replace the actual sample IDs in the `Sample_Name` column if necessary, and to add the sex values to the `Sex` column. Don't change these column names though.

Perform the background correction, dye bias correction and sex prediction:

	qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

We can now summarise the QC analysis of the raw data

	qc.summary <- meffil.qc.summary(qc.objects)

and generate a report:
	
	meffil.qc.report(qc.summary)

This creates a html file, by default called `meffil.qc.report.html` in the current work directory. Should open up in your web browser.

You can remove bad samples prior to performing quantile normalization:

	qc.objects <- meffil.remove.ids(qc.objects, qc.summary$bad.ids)

And now perform quantile normalization:

	qc.quantiles <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

Finally, the `beta` values can be generated, whilst removing CpGs that were found to be dodgy in the QC analysis:

	normalized.beta <- meffil.normalize.samples(qc.quantiles, cpglist.remove=qc.report$bad.cpgs$name)

A summary report of the normalization performance can also be generated:

	normalization.summary <- meffil.normalization.summary(normalized.beta, qc.quantiles)
	meffil.normalization.report(normalization.summary)

By default this will create a file called `meffil.normalization.report.html`.



