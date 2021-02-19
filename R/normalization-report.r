#' Generate report on normalization performance
#'
#' Generate HTML file that summarises the normalization. 
#'
#' @param  normalization.summary Output from \code{meffil.normalization.summary}.
#' @param  output.file Default = "meffil-normalization-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "Illumina methylation data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' @examples \dontrun{
#'
#'}
meffil.normalization.report <- function(
    normalization.summary,
    output.file = "normalization-report.md",
    author = "Analyst",
    study = "Illumina methylation data",
    ...
) {
    msg("Writing report as html file to", output.file)
    report.path <- system.file("reports", package="meffil")
    require(knitr)
    require(Cairo)
    require(gridExtra)
    opts <- opts_chunk$get()
    on.exit(opts_chunk$set(opts))
    opts_chunk$set(warning=FALSE, echo=FALSE, message=FALSE, results="asis", fig.width=6, fig.height=6, dev="CairoPNG")    
    knit.report(file.path(report.path, "normalization-report.rmd"), output.file, ...)
}


#' Perform tests to check normalization performance
#'
#' Creates scree plot of PCs of control probes, tests for association of control probe PCs with batch variables, tests for association of normalized probes with batch variables, creates PCA plots
#' @param  norm.objects Output from \link{meffil.normalize.quantiles}
#' @param  pcs Output from \code{\link{meffil.methylation.pcs}()}
#' applied to the normalized methylation matrix
#' corresponding to \code{norm.objects}.
#' @param  parameters Default = meffil.post.parameters(norm.objects). Report parameters.
#' @param  variables Default = NULL. Data frame of variables to compare to
#' principal components (\code{pcs}).  Must contain \code{length(norm.objects)} rows.
#' Columns that are not factors are ignored.
#' @param  verbose Default = TRUE
#' @export 
#' @return List of tables and graphs.
#' @examples \dontrun{
#'
#'}
meffil.normalization.summary <- function(norm.objects, pcs, parameters = meffil.normalization.parameters(norm.objects), variables=NULL, verbose=TRUE)
{
    stopifnot(sapply(norm.objects, is.normalized.object))
    stopifnot(is.matrix(pcs) && nrow(pcs) == length(norm.objects))
    stopifnot(is.numeric(pcs))
    stopifnot(is.null(variables) || is.data.frame(variables) && nrow(variables) == length(norm.objects))
    
    scree.plot <- meffil.plot.control.scree(norm.objects)
    control.batch <- meffil.plot.control.batch(
        norm.objects, 
        npcs=parameters$control.pcs, 
        variables=parameters$variables,
        additional=variables,
        batch.threshold=parameters$batch.threshold,
        cols=parameters$colours,
        verbose=verbose
	)
    probe.batch <- meffil.plot.probe.batch(
        norm.objects,
        pcs[,intersect(1:ncol(pcs), parameters$batch.pcs),drop=F],
        variables=parameters$variables,
        additional=variables,
        batch.threshold=parameters$batch.threshold,
        cols=parameters$colours,
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
        if (sum(sex == "M") >= 2)
            pcs$male <- meffil.pcs(norm.objects[sex=="M"])
        if (sum(sex == "F") >= 2)
            pcs$female <- meffil.pcs(norm.objects[sex=="F"])
    }
    pc.var <- lapply(pcs, function(pcs) apply(pcs$x, 2, var))
    list(graphs=lapply(names(pcs), function(nom)
             {
                 dat <- data.frame(x=1:length(pc.var[[nom]]), y=pc.var[[nom]])
                 (ggplot(dat, aes(x=x,y=y)) +
                  geom_bar(stat="identity", position="dodge") +
                  labs(x="PC", y="variance") +
                  ggtitle(paste(nom, "samples")))                 
             }),
         tab=pc.var)
}



#' Test for association of control matrix probes with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from the control matrix against each of \code{m} measured batch variables
#'
#' @param  norm.objects From \code{meffil.normalize.quantiles}
#' @param  pcs Which PCs to plot. Default first 10
#' @param  variables. Default = guess.batch.vars(norm.objects). Array spacifying column names in samplesheet to test for association with control matrix PCs.
#' @param additional. Default = NULL. Data frame containing variables to test for association
#' with control matrix PCs. Must have \code{nrow(additional) == length(norm.objects)}.
#' @param  verbose=T Print progress messages?
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.control.batch <- function(norm.objects, npcs=1:10, variables=guess.batch.vars(norm.objects), additional=NULL, batch.threshold=1e-50, cols=NULL, verbose=TRUE)
{
    stopifnot(sapply(norm.objects, is.normalized.object))
    stopifnot(is.null(additional) || is.data.frame(additional) && nrow(additional) == length(norm.objects))

    pcs <- meffil.pcs(norm.objects)$x
    npcs <- intersect(npcs, 1:ncol(pcs))
    stopifnot(length(npcs) > 0)
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
    if (!is.null(additional))
        dat <- cbind(dat, additional)

    msg("Testing associations", verbose=verbose)
    res <- test.pairwise.associations(pcs, dat)
    ret <- plot.pairwise.associations(res, batch.threshold)
    ret$pc.plots <- plot.pcs(pcs, dat, cols)

    colnames(res)[which(colnames(res) == "x")] <- "batch.variable"
    colnames(res)[which(colnames(res) == "l")] <- "batch.value"
    colnames(res)[which(colnames(res) == "y")] <- "principal.component"
    ret$tab <- res
    
    ret
}


#' Test normalized betas for association with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most variable CpG sites (after normalization) against each of \code{m} measured batch variables
#'
#' @param  norm.objects Output from \code{\link{meffil.normalize.quantiles}()}.
#' @param  pcs Output from \code{\link{meffil.methylation.pcs}()}
#' applied to the normalized methylation matrix
#' corresponding to \code{norm.objects}.
#' @param  variables Default = guess.batch.vars(norm). Which variables in sample sheet to test
#' @param additional. Default = NULL. Data frame containing variables to test for association
#' with control matrix PCs. Must have \code{nrow(additional) == length(norm.objects)}.
#' @param  verbose=T Print progress messages?
#' @return List of table of results and graph
#' @examples \dontrun{
#'
#'}
#' @export
meffil.plot.probe.batch <- function(norm.objects, pcs, variables=guess.batch.vars(norm.objects), additional=NULL, batch.threshold=1e-50, cols=NULL, verbose=T)
{
    stopifnot(sapply(norm.objects, is.normalized.object))
    stopifnot(is.matrix(pcs) && nrow(pcs) == length(norm.objects))
    stopifnot(is.numeric(pcs))
    stopifnot(is.null(additional) || is.data.frame(additional) && nrow(additional) == length(norm.objects))

    msg("Extracting batch variables", verbose=verbose)
    variables <- variables[variables %in% names(norm.objects[[1]]$samplesheet)]
    
    dat <- as.data.frame(lapply(variables, function(x) {
        sapply(norm.objects, function(y)
               {
                   y$samplesheet[[x]]
               })
    }), stringsAsFactors=F)
    names(dat) <- variables

    if (!is.null(additional))
        dat <- cbind(dat, additional)
    
    msg("Testing associations", verbose=verbose)
    res <- test.pairwise.associations(pcs, dat)
    ret <- plot.pairwise.associations(res, batch.threshold)
    ret$pc.plots <- plot.pcs(pcs, dat, cols)

    colnames(res)[which(colnames(res) == "x")] <- "batch.variable"
    colnames(res)[which(colnames(res) == "l")] <- "batch.value"
    colnames(res)[which(colnames(res) == "y")] <- "principal.component"
    ret$tab <- res
    ret
}

# Tests associations between 
# between each column of y and each column of x.
# Each column of y must be numeric.
# In cases where a column of x contains factors,
# the y-values for each factor level is compared
# to the set of all y-values that are not outliers.
# Outliers are identified as those below 1.5
# times the inter-quartile range below the first quartile
# or 1.5 times the inter-quartile range above the third quartile.
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
            pval <- NA
            fstat <- NA
            r2 <- NA
            try({
                fit <- summary(lm(y~x))
                fstat <- fit$fstatistic
                if (is.null(fstat)) fstat <- NA
                pval <- pf(fstat["value"], df1=fstat["numdf"], df2=fstat["dendf"], lower.tail=F)#
                r2 <- fit$r.squared
            }, silent=TRUE)
            ret <- data.frame(x=x.name, l=NA, y=y.name, test="F-test",
                              p.value=pval, estimate=fstat["value"], lower=NA, upper=NA, r2=r2)
            if (is.factor(x) && length(levels(x)) > 1) {
                q1 <- quantile(y, probs=0.25)
                q3 <- quantile(y, probs=0.75)
                is.outlier <- y < q1 - 1.5*(q3-q1) | y > q3 + 1.5*(q3-q1)
                x.levels <- names(which(table(x) > 0))
                fits <- sapply(x.levels, function(level) {
                    fit <- NA
                    try({
                        idx <- which(!is.outlier | x == level)
                        fit <- lm(y ~ level, data=data.frame(y=y,level=sign(x==level))[idx,])
                    }, silent=TRUE)
                    fit
                }, simplify=F)
                
                pvals <- sapply(fits, function(fit) {
                    coef <- coefficients(summary(fit))
                    if (is.matrix(coef) && "level" %in% rownames(coef))
                        coef["level","Pr(>|t|)"]
                    else NA
                })

                r2s <- sapply(fits, function(fit) {
                    fit <- summary(fit)
                    if (class(fit) == "summary.lm")
                        fit$r.squared
                    else NA
                })
                                
                confint <- t(sapply(fits, function(fit) {
                    coef <- coefficients(summary(fit))
                    ret <- c(Estimate=NA,lwr=NA,upr=NA)
                    if (is.matrix(coef) && "level" %in% rownames(coef))
                        try(ret <- confint(glht(fit))$confint["level",c("Estimate","lwr","upr")], silent=TRUE)
                    ret
                }))
                                
                colnames(confint) <- c("estimate","lower","upper")
                ret <- rbind(ret,
                             data.frame(x=x.name,
                                        l=x.levels,
                                        y=y.name,
                                        test="t-test",
                                        p.value=pvals,
                                        confint,
                                        r2=r2s))
            }
            ret
        }))
    }))
    rownames(ret) <- NULL
    ret
}




plot.pairwise.associations <- function(res, batch.threshold) {
    res <- res[which(res$test == "F-test" | res$test == "t-test" & !is.na(res$estimate)),]
    
    cp <- sapply(unique(res$y), function(y.value) {
        idx <- which(res$y == y.value & res$test == "t-test")  
        col <- with(res[idx,], c("black","red")[sign(p.value < batch.threshold) + 1])
        (ggplot(res[idx,], aes(x=paste(x, l, sep="."),
                               y=estimate,
                               ymin=lower,
                               ymax=upper)) +
         geom_errorbar(width=0.5, colour=col) +
         geom_point(colour=col) +
         geom_hline(yintercept=0, colour="blue", linetype="dashed") +
         ylab("coefficient") +
         xlab("") +
         ggtitle(y.value) +
         theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
         coord_flip())
    }, simplify=F)

    res <- res[which(res$test == "F-test"),]

    fp <- sapply(unique(res$x), function(x.value) {
        idx <- which(res$x == x.value)
        (ggplot(res[idx,], aes(x=y, y=-log10(p.value))) +
         geom_point() +
         geom_hline(yintercept=-log10(0.05), colour="blue", linetype="dashed") +
         labs(x="", y="-log10 p") +
         ggtitle(x.value) +
         theme_bw())
    }, simplify=F)
    
    return(list(fplots=fp, cplots=cp))
}


plot.pcs <- function(pcs, dat, cols=NULL) {
    stopifnot(nrow(pcs) == nrow(dat))

    if (is.null(cols)) {
        ## http://colorbrewer2.org/
        ## 12 colours, qualitative
        cols <- c("#a6cee3", "#1f78b4", "#b2df8a",
                  "#33a02c","#fb9a99", "#e31a1c",
                  "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a","#ffff99", "#b15928")
    }
    
    if (ncol(pcs) >= 3) {
        too.many.levels <- sapply(1:ncol(dat), function(i) length(unique(dat[,i])) > 10)
        if (any(!too.many.levels)) {            
            pc.vars <- do.call(rbind, lapply(which(!too.many.levels), function(i) {                
                rbind(data.frame(desc="pc1vpc2", pc.x=pcs[,1], pc.y=pcs[,2],
                                 variable=colnames(dat)[i],
                                 values=paste(colnames(dat)[i], dat[,i], sep="."),
                                 stringsAsFactors=F),
                      data.frame(desc="pc1vpc3", pc.x=pcs[,1], pc.y=pcs[,3],
                                 variable=colnames(dat)[i],
                                 values=paste(colnames(dat)[i], dat[,i], sep="."),
                                 stringsAsFactors=F),
                      data.frame(desc="pc2vpc3", pc.x=pcs[,2], pc.y=pcs[,3],
                                 variable=colnames(dat)[i],
                                 values=paste(colnames(dat)[i], dat[,i], sep="."),
                                 stringsAsFactors=F))
            }))

            return(sapply(unique(pc.vars$variable), function(variable) {
                idx <- which(pc.vars$variable == variable)
                n.values <- length(unique(pc.vars$values[idx]))
                if (n.values > length(cols))
                    cols <- rep(cols, length.out=n.values)
                
                (ggplot(pc.vars[idx,], aes(x=pc.x, y=pc.y,colour=as.factor(values))) +
                 geom_point() +
                 scale_colour_manual(name="Batch", values=cols) +
                 labs(y="pc",x="pc") +
                 facet_grid(variable ~ desc) +
                 theme_bw())
            },simplify=F))            
        }
    }
    return(NULL)
}


#' Guess which columns in sample sheet are batch variables
#'
#' @param norm.objects Output from \link{meffil.normalize.quantiles}
#' @export
#' @return Array of variable names
#' @examples \dontrun{
#'
#'}
guess.batch.vars <- function(norm.objects) {
    nom <- names(norm.objects[[1]]$samplesheet)
    nom <- nom[!nom %in% c("Sample_Name", "Sex", "Basename")]
    return(nom)
}


#' Specify parameters for testing normalization
#'
#' @param  norm.objects Output from \link{meffil.normalize.quantiles}
#' @param  variables Default = guess.batch.vars(norm). Which variables in sample sheet to test
#' @param  control.pcs Default = 1:10. Number of control PCs to test against batch variables
#' @param  colours Colours to use for scatterplots.
#' @export
#' @return List of parameters
#' @examples \dontrun{
#'
#'}
meffil.normalization.parameters <- function(norm.objects,
                                            variables = guess.batch.vars(norm.objects),
                                            control.pcs = 1:10,
                                            batch.pcs = 1:10,
                                            batch.threshold=1e-50,
                                            colours=NULL) {                                           
    stopifnot(sapply(norm.objects, is.normalized.object))
    
    parameters <- list(
        variables = variables,
        control.pcs = control.pcs,
        batch.pcs = batch.pcs,
        batch.threshold=batch.threshold,
        colours=colours
	)
    return(parameters)
}


