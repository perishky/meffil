library(ggplot2)
library(reshape2)

#' Plot scree plot of control matrix
#'
#' @param  qc.quantiles From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}

meffil.plot.control.scree <- function(qc.quantiles)
{
	nom <- names(qc.quantiles$pca)
	par(mfrow=c(length(nom), 1))
	for(i in 1:length(nom))
	{
		dat <- qc.quantiles$pca[[i]]$sdev^2
		names(dat) <- 1:length(dat)
		barplot(dat, main=nom[i], ylab="Variance", xlab="PC")
	}
}


#' Test for association of control matrix probes with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most control matrix against each of \code{m} measured batch variables
#'
#' @param  qc.quantiles From \code{meffil.normalize.objects}
#' @param  pcs Which PCs to plot. Default first 10
#' @param  variables. Default = guess.batch.vars(qc.quantiles). Array spacifying column names in samplesheet to test for association with control matrix PCs
#' @param  verbose=T Print progress messages?
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.control.batch <- function(qc.quantiles, npcs=1:10, variables=guess.batch.vars(qc.quantiles), verbose=TRUE)
{
	stopifnot(all(npcs %in% 1:ncol(qc.quantiles$pca$all$x)))
	pcs <- qc.quantiles$pca$all$x[,npcs]
	objects <- qc.quantiles$samples
	variables <- variables[variables %in% names(objects[[1]]$samplesheet)]

	msg("Extracting batch variables", verbose=verbose)
	dat <- as.data.frame(lapply(variables, function(x) {
		sapply(objects, function(y)
		{
			y$samplesheet[[x]]
		})
	}))
	names(dat) <- variables

	msg("Testing associations", verbose=verbose)
	res <- matrix(0, ncol(dat), ncol(pcs))
	for(i in 1:ncol(dat))
	{
		for(j in 1:ncol(pcs))
		{
			res[i,j] <- coefficients(summary(lm(pcs[,j] ~ dat[,i])))[2,4]
		}
	}
	colnames(res) <- paste("PC", 1:ncol(pcs), sep="")
	res <- data.frame(v=variables, res)
	res <- reshape2::melt(res, id.vars="v", measure.vars=paste("PC", 1:ncol(pcs), sep=""))
	p1 <- ggplot(res, aes(x=variable, y=-log10(value))) +
	geom_point() +
	geom_hline(yintercept=-log10(0.05), linetype="dotted") +
	facet_grid(v ~ .) +
	labs(y="-log10 p", x="PCs") +
	theme_bw()
	return(list(tab=res, graph=p1))
}


#' Test normalized betas for association with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most variable CpG sites against each of \code{m} measured batch variables
#'
#' @param  normalized.beta Output from \code{meffil.normalize.samples}
#' @param  qc.quantiles Output from \code{meffil.normalize.quantiles}
#' @param  npcs Default = 1:10. Which CpG PCs to test
#' @param  variables Default = guess.batch.vars(qc.quantiles). Which variables in sample sheet to test
#' @param  probe.range Default = 1:1000. How many probes to be used in calculating PCs
#' @param  verbose=T Print progress messages?
#' @export
#' @return List of table of results and graph
#' @examples \dontrun{
#'
#'}
meffil.plot.probe.batch <- function(normalized.beta, qc.quantiles, npcs=1:10, variables=guess.batch.vars(qc.quantiles), probe.range=1000, verbose=T)
{
	stopifnot(all(npcs %in% 1:ncol(normalized.beta)) & all(npcs %in% 1:probe.range))

	msg("Calculating variances", verbose=verbose)
	vars <- apply(normalized.beta, 1, var)
	varids <- order(vars, decreasing=TRUE)[1:probe.range]

	msg("Calculating beta PCs", verbose=verbose)
	pcs <- prcomp(t(normalized.beta[varids,]))$x[,npcs]

	msg("Extracting batch variables", verbose=verbose)
	objects <- qc.quantiles$samples
	variables <- variables[variables %in% names(objects[[1]]$samplesheet)]

	dat <- as.data.frame(lapply(variables, function(x) {
		sapply(objects, function(y)
		{
			y$samplesheet[[x]]
		})
	}))
	names(dat) <- variables

	msg("Testing associations", verbose=verbose)
	res <- matrix(0, ncol(dat), ncol(pcs))
	for(i in 1:ncol(dat))
	{
		for(j in 1:ncol(pcs))
		{
			res[i,j] <- coefficients(summary(lm(pcs[,j] ~ dat[,i])))[2,4]
		}
	}
	colnames(res) <- paste("PC", 1:ncol(pcs), sep="")
	res <- data.frame(v=variables, res)
	res <- reshape2::melt(res, id.vars="v", measure.vars=paste("PC", 1:ncol(pcs), sep=""))
	p1 <- ggplot(res, aes(x=variable, y=-log10(value))) +
	geom_point() +
	geom_hline(yintercept=-log10(0.05), linetype="dotted") +
	facet_grid(v ~ .) +
	labs(y="-log10 p", x="PCs") +
	theme_bw()
	return(list(tab=res, graph=p1))

}


#' Guess which columns in sample sheet are batch variables
#'
#' @param qc.quantiles Output from \link{meffil.normalize.quantiles}
#' @export
#' @return Array of variable names
#' @examples \dontrun{
#'
#'}
guess.batch.vars <- function(qc.quantiles)
{
	nom <- names(qc.quantiles$samples[[1]]$samplesheet)
	nom <- nom[!nom %in% c("Sample_Name", "Sex", "Basename")]
	return(nom)
}


#' Specify parameters for testing normalization
#'
#' @param  qc.quantiles Output from \link{meffil.normalize.quantiles}
#' @param  variables Default = guess.batch.vars(qc.quantiles). Which variables in sample sheet to test
#' @param  control.pcs Default = 1:10. Number of control PCs to test against batch variables
#' @param  probe.pcs Default = 1:10. Number of probe PCs to test against batch variables
#' @param  probe.range Default = 1:1000. Number of probes to use to calculate PCs for actual CpGs
#' @export
#' @return List of parameters
#' @examples \dontrun{
#'
#'}
meffil.normalization.parameters <- function(qc.quantiles, variables = guess.batch.vars(qc.quantiles), control.pcs = 1:10, probe.pcs = 1:10, probe.range = 1000)
{
	parameters <- list(
		variables = variables,
		control.pcs = control.pcs,
		probe.pcs = probe.pcs,
		probe.range = probe.range
	)
	return(parameters)
}


#' Perform tests to check normalization performance
#'
#' Creates scree plot of PCs of control probes, tests for association of control probe PCs with batch variables, tests for association of normalized probes with batch variables
#' @param  normalized.beta Output from \code{meffil.normalize.samples}
#' @param  qc.quantiles Output from \link{meffil.normalize.quantiles}
#' @param  parameters Default = meffil.normalization.parameters(qc.quantiles) List of parameters.
#' @param  verbose Default = TRUE
#' @export 
#' @return List of tables and graphs
#' @examples \dontrun{
#'
#'}
meffil.normalization.summary <- function(normalized.beta, qc.quantiles, parameters = meffil.normalization.parameters(qc.quantiles), verbose=TRUE)
{
	scree.plot <- meffil.plot.control.scree(qc.quantiles)
	control.batch <- meffil.plot.control.batch(
		qc.quantiles, 
		npcs=parameters$control.pcs, 
		variables=parameters$variables, 
		verbose=verbose
	)
	probe.batch <- meffil.plot.probe.batch(
		normalized.beta,
		qc.quantiles,
		npcs=parameters$probe.pcs,
		variables=parameters$variables,
		probe.range=parameters$probe.range
	)
	return(list(
		scree.plot = scree.plot,
		control.batch = control.batch,
		probe.batch = probe.batch,
		parameters = parameters
	))
}


#' Generate report on normalization performance
#'
#' Generate HTML file that summarises the normalization. 
#'
#' @param  normalization.summary Output from \code{meffil.normalization.summary}.
#' @param  output.file Default = "meffil.normalization.report.html"
#' If specified then a html report will be generated summarising the normalization.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  studyname Default = "IlluminaHuman450 data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{rmarkdown::render}}
#' @export
#' @return NULL
#' @examples \dontrun{
#'
#'}
meffil.normalization.report <- function(
    normalization.summary,
    output.file = "meffil.normalization.report.html",
    author = "Analyst",
    studyname = "IlluminaHuman450 data",
    ...
) {
    cat("Writing report as html file to", output.file, "\n")
    save(normalization.summary,
        author,
        studyname,
        file = file.path(tempdir(), "meffil.normalization.report.rdata")
    )
    output.dir <- ifelse(dirname(output.file) == ".", getwd(), dirname(output.file))
    rmarkdown::render(system.file("reports", "meffil.normalization.report.rmd", package="meffil"), output_file=basename(output.file), output_dir=output.dir, ...)
}
