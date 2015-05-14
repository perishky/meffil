
#' Generate report on normalization performance
#'
#' Generate HTML file that summarises the normalization. 
#'
#' @param  normalization.summary Output from \code{meffil.normalization.summary}.
#' @param  output.file Default = "meffil-normalization-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "IlluminaHuman450 data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' @examples \dontrun{
#'
#'}
meffil.normalization.report <- function(
    normalization.summary,
    output.file = "meffil-normalization-report.md",
    author = "Analyst",
    study = "IlluminaHuman450 data",
    ...
) {
    msg("Writing report as html file to", output.file)
    path <- system.file("reports", package="meffil")
    knit.report(file.path(path, "meffil-normalization-report.rmd"), output.file, ...)
}


#' Perform tests to check normalization performance
#'
#' Creates scree plot of PCs of control probes, tests for association of control probe PCs with batch variables, tests for association of normalized probes with batch variables
#' @param  normalized.beta Output from \code{meffil.normalize.samples}
#' @param  norm Output from \link{meffil.normalize.quantiles}
#' @param  parameters Default = meffil.post.parameters(norm) List of parameters.
#' @param  verbose Default = TRUE
#' @export 
#' @return List of tables and graphs
#' @examples \dontrun{
#'
#'}
meffil.normalization.summary <- function(normalized.beta, norm.objects, parameters = meffil.normalization.parameters(norm.objects), verbose=TRUE)
{
    stopifnot(sapply(norm.objects, is.normalization.object))
    
    scree.plot <- meffil.plot.control.scree(norm.objects)
    control.batch <- meffil.plot.control.batch(
        norm.objects, 
        npcs=parameters$control.pcs, 
        variables=parameters$variables, 
        verbose=verbose
	)
    probe.batch <- meffil.plot.probe.batch(
        normalized.beta,
        norm.objects,
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



#' Plot scree plot of control matrix
#'
#' @param  norm.objects From \code{meffil.normalize.quantiles}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}

meffil.plot.control.scree <- function(norm.objects)
{
    stopifnot(sapply(norm.objects, is.normalization.object))

    sex <- sapply(norm.objects, function(object) object$predicted.sex)
    pcs <- list(all=meffil.pcs(norm.objects))
    if (max(table(sex)) < length(sex)) {
        pcs$male <- meffil.pcs(norm.objects[sex=="M"])
        pcs$female <- meffil.pcs(norm.objects[sex=="F"])
    }
    pc.var <- lapply(pcs, function(pcs) apply(pcs$x, 2, var))
    list(graphs=lapply(names(pcs), function(nom)
             {
                 qplot(x=1:length(pc.var[[nom]]), y=pc.var[[nom]],
                       geom="bar", stat="identity", position="dodge") +
                       labs(x="PC", y="variance") +
                       ggtitle(paste(nom, "samples"))                           
             }),
         tab=pc.var)
}


#' Test for association of control matrix probes with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most control matrix against each of \code{m} measured batch variables
#'
#' @param  norm.objects From \code{meffil.normalize.quantiles}
#' @param  pcs Which PCs to plot. Default first 10
#' @param  variables. Default = guess.batch.vars(norm.objects). Array spacifying column names in samplesheet to test for association with control matrix PCs
#' @param  verbose=T Print progress messages?
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.control.batch <- function(norm.objects, npcs=1:10, variables=guess.batch.vars(norm.objects), verbose=TRUE)
{
    stopifnot(sapply(norm.objects, is.normalization.object))

    pcs <- meffil.pcs(norm.objects)$x
    stopifnot(all(npcs %in% 1:ncol(pcs)))
    pcs <- pcs[,npcs,drop=F]
    
    variables <- variables[variables %in% names(norm.objects[[1]]$samplesheet)]
    
    msg("Extracting batch variables", verbose=verbose)
    dat <- as.data.frame(lapply(variables, function(x) {
        sapply(norm.objects, function(y)
               {
                   y$samplesheet[[x]]
               })
    }), stringsAsFactors=F)
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
    res <- data.frame(v=variables, res, stringsAsFactors=F)
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
#' @param  norm.objects Output from \code{meffil.normalize.quantiles}
#' @param  npcs Default = 1:10. Which CpG PCs to test
#' @param  variables Default = guess.batch.vars(norm). Which variables in sample sheet to test
#' @param  probe.range Default = 1:1000. How many probes to be used in calculating PCs
#' @param  verbose=T Print progress messages?
#' @export
#' @return List of table of results and graph
#' @examples \dontrun{
#'
#'}
meffil.plot.probe.batch <- function(normalized.beta, norm.objects, npcs=1:10, variables=guess.batch.vars(norm.objects), probe.range=1000, verbose=T)
{
    stopifnot(sapply(norm.objects, is.normalization.object))
    
    stopifnot(all(npcs %in% 1:ncol(normalized.beta)) & all(npcs %in% 1:probe.range))
    
    msg("Calculating variances", verbose=verbose)
    
    if (class(normalized.beta) == "big.matrix") {
        subset.idx <- sample(1:nrow(normalized.beta), size=floor(nrow(normalized.beta))*0.1)
        normalized.beta <- normalized.beta[subset.idx,]
    }
    vars <- matrixStats::rowVars(normalized.beta, na.rm=T)
    varids <- order(vars, decreasing=TRUE)[1:probe.range]
    
    msg("Calculating beta PCs", verbose=verbose)
    pcs <- prcomp(t(normalized.beta[varids,]))$x[,npcs,drop=F]
    
    msg("Extracting batch variables", verbose=verbose)
    variables <- variables[variables %in% names(norm.objects[[1]]$samplesheet)]
    
    dat <- as.data.frame(lapply(variables, function(x) {
        sapply(norm.objects, function(y)
               {
                   y$samplesheet[[x]]
               })
    }), stringsAsFactors=F)
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
    res <- data.frame(v=variables, res, stringsAsFactors=F)
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
#' @param norm.objects Output from \link{meffil.normalize.quantiles}
#' @export
#' @return Array of variable names
#' @examples \dontrun{
#'
#'}
guess.batch.vars <- function(norm.objects)
{
	nom <- names(norm.objects[[1]]$samplesheet)
	nom <- nom[!nom %in% c("Sample_Name", "Sex", "Basename")]
	return(nom)
}


#' Specify parameters for testing normalization
#'
#' @param  norm.objects Output from \link{meffil.normalize.quantiles}
#' @param  variables Default = guess.batch.vars(norm). Which variables in sample sheet to test
#' @param  control.pcs Default = 1:10. Number of control PCs to test against batch variables
#' @param  probe.pcs Default = 1:10. Number of probe PCs to test against batch variables
#' @param  probe.range Default = 1:1000. Number of probes to use to calculate PCs for actual CpGs
#' @export
#' @return List of parameters
#' @examples \dontrun{
#'
#'}
meffil.normalization.parameters <- function(norm.objects, variables = guess.batch.vars(norm.objects), control.pcs = 1:10, probe.pcs = 1:10, probe.range = 1000)
{
    stopifnot(sapply(norm.objects, is.normalization.object))
    
    parameters <- list(
        variables = variables,
        control.pcs = control.pcs,
        probe.pcs = probe.pcs,
        probe.range = probe.range
	)
    return(parameters)
}


