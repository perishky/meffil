# meffil
Efficient algorithms for analyzing DNA methylation data
generated using Infinium HumanMethylation450 or MethylationEPIC BeadChips:

* Functional normalization for large datasets using parallelization.
* Normalization of datasets with mixed Infinium HumanMethylation450 and MethylationEPIC BeadChips.
* Inclusion of user-defined fixed and random effects in functional normalization procedure.
* Cell count estimation using predefined and user-defined reference datasets.
* Use of predefined and user-defined microarray probe annotations.
* Epigenome-wide association studies (using data from any normalization pipeline).
* Copy number estimation.
* Report generation summarizing all steps.

The `minfi` version consists of a single function call that requires
access to all data. The `meffil` version splits the method up into
several functions in order to allow parallelization and to reduce the
amount of data that needs to be loaded.

## One-step normalization

	library(meffil)
	options(mc.cores=6)

	# Generate samplesheet
	samplesheet <- meffil.create.samplesheet(path_to_idat_files)

	# Or read in samplesheet
	samplesheet <- meffil.read.samplesheet(path_to_idat_files)

    beta <- meffil.normalize.dataset(samplesheet, qc.file="qc/report.html", author="Analyst", study="Illumina450", number.pcs=10)

## Step-by-step normalization

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
	meffil.qc.report(qc.summary, output.file="qc/report.html")

	# Remove outlier samples if necessary
	qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

    # Plot residuals remaining after fitting control matrix to decide on the number PCs
	# to include in the normalization below.
	print(meffil.plot.pc.fit(qc.objects)$plot)

	# Perform quantile normalization
	norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

	# Generate normalized probe values
	norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)

	# Generate normalization report
	norm.summary <- meffil.normalization.summary(norm.beta, norm.objects)
	meffil.normalization.report(norm.summary, output.file="normalization/report.html")

## More info about normalization

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

These data don't actually have a samplesheet, so we can generate one from the `.idat` files:

	samplesheet <- meffil.create.samplesheet(path)

The function creates the following necessary columns:

- Sample_Name
- Sex (possible values "M"/"F"/NA)
- Basename

And it also tries to parse the basenames to guess if the Sentrix plate
and positions are present. At this point it is worthwhile to manually
modify the `samplesheet` data.frame to replace the actual sample IDs
in the `Sample_Name` column if necessary, and to add the sex values to
the `Sex` column. Don't change these column names though.

Perform the background correction, dye bias correction, sex prediction and cell count estimates:

	qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE)

A list of available cell type references can be obtained as follows:

	meffil.list.cell.type.references()

New references can be created from a dataset using meffil.create.cell.type.reference().

Obtain the matrix of genotypes for comparison with those measured on the microarray.
If such a matrix is available (rows = SNPs, columns = samples), then the following steps
can be omitted.  Otherwise, it is possible to obtain the matrix from a PLINK
dataset as follows:

    annotation <- qc.objects[[1]]$annotation
    writeLines(meffil.snp.names(annotation), con="snp-names.txt")
    command shell > plink --bfile dataset --extract snp-names.txt --recodeA --out genotypes.raw --noweb
    filenames <- "genotypes.raw"
    genotypes <- meffil.extract.genotypes(filenames)

We can now summarise the QC analysis of the raw data

	qc.summary <- meffil.qc.summary(qc.objects, genotypes=genotypes)

and generate a report:
	
	meffil.qc.report(qc.summary, output.file="qc/report.html")

This creates the file "qc/report.html" in the current work directory.
Should open up in your web browser.

You can remove bad samples prior to performing quantile normalization:

	qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

Next we determine the number of principal components of the control matrix
to include in the quantile normalization.
The following function plots the quantile residuals remaining
after fitting different numbers of control matrix principal components.

    print(meffil.plot.pc.fit(qc.objects)$plot)

And now remove control probe variance from the sample quantiles:

	norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

Additional fixed and random effects can be included in the normalization
by providing their corresponding column names in the samplesheet.
For example, slide effects can be included as follows:

    norm.objects <- meffil.normalize.quantiles(qc.objects, random.effects="Slide", number.pcs=10)

Note however that including random effects will greatly increase running time.

Finally, the `beta` values can be generated, whilst removing CpGs
that were found to be dodgy in the QC analysis:

	norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)

A summary report of the normalization performance can also be generated:

	norm.summary <- meffil.normalization.summary(norm.beta, norm.objects)
	meffil.normalization.report(norm.summary, output.file="normalization/report.html")

## Epigenome-wide association study (EWAS)

Prepare a variable of interest and a set of covariates.

    variable <- ... ## variable of interest, one value per sample
    covariates <- ... ## data.frame of covariates to include (rows = samples, columns = covariates)

Add cell count estimates to the set  of covariates.

    counts <- t(meffil.cell.count.estimates(norm.objects))
	covariates <- cbind(covariates, counts)

Run the EWAS.

    ewas.ret <- meffil.ewas(norm.beta, variable=variable, covariates=covariates)

Generate a report for the EWAS, including a table describing
all variables, a table describing the relationships between the variable of interest
and all covariates, QQ plots, Manhattan plots, a spreadsheet listing all significantly associated CpG sites,
plots of the most strongly associated sites,
and plots of selected candidate CpG sites.

    ewas.parameters <- meffil.ewas.parameters(sig.threshold=1e-20, max.plots=5)
    candidate.sites <- c("cg04946709","cg06710937","cg12177922","cg15817705","cg20299935","cg21784396")
    ewas.summary <- meffil.ewas.summary(ewas.ret,
                                        norm.beta,
                                        selected.cpg.sites=candidate.sites,
                                        parameters=ewas.parameter)								
    meffil.ewas.report(ewas.summary, output.file="ewas/report.html")




