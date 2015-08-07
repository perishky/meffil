
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
#' Creates scree plot of PCs of control probes, tests for association of control probe PCs with batch variables, tests for association of normalized probes with batch variables, creates PCA plots
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
    stopifnot(sapply(norm.objects, is.normalized.object))
    
    scree.plot <- meffil.plot.control.scree(norm.objects)
    control.batch <- meffil.plot.control.batch(
        norm.objects, 
        npcs=parameters$control.pcs, 
        variables=parameters$variables,
        batch.threshold=parameters$batch.threshold,
        verbose=verbose
	)
    probe.batch <- meffil.plot.probe.batch(
        normalized.beta,
        norm.objects,
        npcs=parameters$probe.pcs,
        variables=parameters$variables,
        probe.range=parameters$probe.range,
        batch.threshold=parameters$batch.threshold,
        verbose=verbose
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
    stopifnot(sapply(norm.objects, is.normalized.object))

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
#' Performs association of each of \code{n} PCs calculated from the control matrix against each of \code{m} measured batch variables
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
meffil.plot.control.batch <- function(norm.objects, npcs=1:10, variables=guess.batch.vars(norm.objects), batch.threshold=1e-50, verbose=TRUE)
{
    stopifnot(sapply(norm.objects, is.normalized.object))

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
    res <- test.pairwise.associations(pcs, dat)
    res <- res[which(res$test == "F-test" | res$p.value < batch.threshold),]
    res <- res[order(res$test, res$p.value, decreasing=F),]
    colnames(res)[which(colnames(res) == "x")] <- "batch.variable"
    colnames(res)[which(colnames(res) == "y")] <- "PC"

    msg("Plotting PCs", verbose=verbose)
    scores<-data.frame()
    for (i in 1:ncol(dat)){
    d <- data.frame(PC="PC1vsPC2",v=variables[i],var=dat[,i], PC.scores.x=pcs[,1],PC.scores.y=pcs[,2])
    d2<- data.frame(PC="PC1vsPC3",v=variables[i],var=dat[,i], PC.scores.x=pcs[,1],PC.scores.y=pcs[,3])
    d3<-data.frame(PC="PC2vsPC3",v=variables[i],var=dat[,i], PC.scores.x=pcs[,2],PC.scores.y=pcs[,3])
    d<-rbind(d,d2)
    d<-rbind(d,d3)
    scores <-rbind(scores,d)
    }

    pcaplot <- ggplot(scores,aes(x=PC.scores.x, y=PC.scores.y,colour=factor(scores$var))) +
    geom_point() +
    scale_colour_discrete(name = "Batch") +
    labs(y="PC.scores",x="PC.scores") +
    facet_grid(v ~ PC) +
    theme_bw()

    p1 <- ggplot(res[which(res$test == "F-test"),], aes(x=PC, y=-log10(p.value))) +
	geom_point() +
        geom_hline(yintercept=-log10(0.05), linetype="dotted") +
        facet_grid(batch.variable ~ .) +
        labs(y="-log10 p", x="PCs") +
        theme_bw()
    
    return(list(tab=res, graph=p1,pcaplot=pcaplot))
}


#' Test normalized betas for association with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most variable CpG sites (after normalization) against each of \code{m} measured batch variables
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
meffil.plot.probe.batch <- function(normalized.beta, norm.objects, npcs=1:10, variables=guess.batch.vars(norm.objects), batch.threshold=1e-50, probe.range=5000, verbose=T)
{
    stopifnot(sapply(norm.objects, is.normalized.object))
    
    stopifnot(all(npcs %in% 1:ncol(normalized.beta)) & all(npcs %in% 1:probe.range))
    
    msg("Calculating variances", verbose=verbose)
    
    if (class(normalized.beta) == "big.matrix") {
        subset.idx <- sample(1:nrow(normalized.beta), size=floor(nrow(normalized.beta))*0.1)
        normalized.beta <- normalized.beta[subset.idx,]
    }
    autosomal.idx <- which(rownames(normalized.beta) %in% meffil.get.autosomal.sites())
    vars <- matrixStats::rowVars(normalized.beta[autosomal.idx,], na.rm=T)
    varids <- order(vars, decreasing=TRUE)[1:probe.range]
    
    msg("Calculating beta PCs", verbose=verbose)
    pcs <- prcomp(t(impute.matrix(normalized.beta[varids,], margin=1)))$x[,npcs,drop=F]
    
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
    res <- test.pairwise.associations(pcs, dat)
    res <- res[which(res$test == "F-test" | res$p.value < batch.threshold),]
    res <- res[order(res$test, res$p.value, decreasing=F),]
    colnames(res)[which(colnames(res) == "x")] <- "batch.variable"
    colnames(res)[which(colnames(res) == "y")] <- "PC"

   msg("Plotting PCs", verbose=verbose)
    scores<-data.frame()
    for (i in 1:ncol(dat)){
    d <- data.frame(PC="PC1vsPC2",v=variables[i],var=dat[,i], PC.scores.x=pcs[,1],PC.scores.y=pcs[,2])
    d2<- data.frame(PC="PC1vsPC3",v=variables[i],var=dat[,i], PC.scores.x=pcs[,1],PC.scores.y=pcs[,3])
    d3<-data.frame(PC="PC2vsPC3",v=variables[i],var=dat[,i], PC.scores.x=pcs[,2],PC.scores.y=pcs[,3])
    d<-rbind(d,d2)
    d<-rbind(d,d3)
    scores <-rbind(scores,d)
    }

    pcaplot <- ggplot(scores,aes(x=PC.scores.x, y=PC.scores.y,colour=factor(scores$var))) +
    geom_point() +
    scale_colour_discrete(name = "Batch") +
    labs(y="PC.scores",x="PC.scores") +
    facet_grid(v ~ PC) +
    theme_bw()
    return(list(tab=res, graph=pcaplot))

    p1 <- ggplot(res[which(res$test == "F-test"),], aes(x=PC, y=-log10(p.value))) +
	geom_point() +
	geom_hline(yintercept=-log10(0.05), linetype="dotted") +
	facet_grid(batch.variable ~ .) +
	labs(y="-log10 p", x="PCs") +
	theme_bw()
    return(list(tab=res, graph=p1))
}

#' Tests associations between 
#' between each column of y and each column of x.
#' Each column of y must be numeric.
#' In cases where a column of x contains factors,
#' the y-values for each factor level is compared
#' to the set of all y-values that are not outliers.
#' Outliers are identified as those below 1.5
#' times the inter-quartile range below the first quartile
#' or 1.5 times the inter-quartile range above the third quartile.
test.pairwise.associations <- function(y,x) {
    stopifnot(nrow(x) == nrow(y))
    ret <- do.call(rbind, lapply(1:ncol(y), function(i) {
        y.name <- colnames(y)[i]
        y <- y[,i]
        stopifnot(is.numeric(y))
        do.call(rbind, lapply(1:ncol(x), function(j) {
            x.name <- colnames(x)[j]
            x <- x[,j]
            if (!is.numeric(x)) x <- as.factor(x)            
            pval <- 1
            try({
                fstat <- summary(lm(y ~ x))$fstatistic
                pval <- pf(fstat["value"], df1=fstat["numdf"], df2=fstat["dendf"], lower.tail=F)#
            }, silent=TRUE)

            ret <- data.frame(x=x.name, y=y.name, test="F-test", p.value=pval)
            if (is.factor(x) && length(levels(x)) > 2) {
                q1 <- quantile(y, probs=0.25)
                q3 <- quantile(y, probs=0.75)
                is.outlier <- y < q1 - 1.5*(q3-q1) | y > q3 + 1.5*(q3-q1)
                pvals <- sapply(levels(x), function(level) {
                    pval <- 1
                    try({
                        idx <- which(!is.outlier | x == level)
                        fit <- lm(y ~ l, data=data.frame(y=y,l=sign(x==level))[idx,])
                        pval <- coefficients(summary(fit))["l","Pr(>|t|)"]
                    }, silent=TRUE)
                    pval
                })
                ret <- rbind(ret,
                             data.frame(x=paste(x.name, levels(x), sep="."),
                                        y=y.name,
                                        test="t-test",
                                        p.value=pvals))
            }                
            ret            
        }))
    }))
    rownames(ret) <- NULL
    ret
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
meffil.normalization.parameters <- function(norm.objects, variables = guess.batch.vars(norm.objects), control.pcs = 1:10, probe.pcs = 1:10, probe.range = 5000, batch.threshold=1e-50)
{
    stopifnot(sapply(norm.objects, is.normalized.object))
    
    parameters <- list(
        variables = variables,
        control.pcs = control.pcs,
        probe.pcs = probe.pcs,
        probe.range = probe.range,
        batch.threshold=batch.threshold
	)
    return(parameters)
}


