meffil.normalization.report.from.betas <- function(
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
    knit.report(file.path(report.path, "normalization-report.from.betas.rmd"), output.file, ...)
}


#' Perform tests to check normalization performance
#'
#' @param  pcs Output from \code{\link{meffil.methylation.pcs}()}
#' applied to the normalized methylation matrix
#' corresponding to \code{norm.objects}.
#' @param  parameters Default = meffil.normalization.parameters.from.betas(). Report parameters.
#' @param  samplesheet Default = NULL. Data frame of variables to compare to
#' principal components (\code{pcs}).  Must contain \code{nrow(pcs) == nrow(samplesheet)} rows.
#' Columns that are not factors are ignored.
#' @param  variables Which variables in sample sheet to test
#' @param  verbose Default = TRUE
#' @export 
#' @return List of tables and graphs.
#' @examples \dontrun{
#'
#'}

meffil.normalization.summary.from.betas <- function(pcs, parameters = meffil.normalization.parameters.from.betas(), samplesheet=samplesheet,variables=variables,verbose=TRUE)
{
    
    stopifnot(is.matrix(pcs) && nrow(pcs) == nrow(samplesheet))
    stopifnot(is.numeric(pcs))
    stopifnot(samplesheet$Sample_Name==names(pcs))
  
    probe.batch <- meffil.plot.probe.batch.from.betas(
        samplesheet,
        pcs[,intersect(1:ncol(pcs), parameters$batch.pcs),drop=F],
        variables=variables,
        batch.threshold=parameters$batch.threshold,
        cols=parameters$colours,
        verbose=verbose
    )
    return(list(
        probe.batch = probe.batch,
        parameters = parameters
    ))
}

#' Test normalized betas for association with known batch variables
#'
#' Performs association of each of \code{n} PCs calculated from most variable CpG sites (after normalization) against each of \code{m} measured batch variables
#'
#' @param  norm.objects Output from \code{\link{meffil.normalize.quantiles}()}.
#' @param  pcs Output from \code{\link{meffil.methylation.pcs}()}
#' applied to the normalized methylation matrix
#' corresponding to \code{norm.objects}.
#' @param  variables Which variables in sample sheet to test
#' @param samplesheet. Data frame containing variables to test for association
#' with control matrix PCs. Must have \code{nrow(pcs) == nrow(samplesheet)}.
#' @param  verbose=T Print progress messages?
#' @return List of table of results and graph
#' @examples \dontrun{
#'
#'}
#' @export
meffil.plot.probe.batch.from.betas <- function(samplesheet,variables,pcs, batch.threshold=batch.threshold, cols=NULL, verbose=T)
{
    stopifnot(is.matrix(pcs) && nrow(pcs) == nrow(samplesheet))
    stopifnot(is.numeric(pcs))
    stopifnot(samplesheet$Sample_Name==names(pcs))
    
    msg("Extracting batch variables", verbose=verbose)
    
    msg("Testing associations", verbose=verbose)
    dat<-samplesheet[,which(names(samplesheet)%in%variables)]
    res <- test.pairwise.associations(pcs, dat)

    ret <- plot.pairwise.associations(res, batch.threshold)
    ret$pc.plots <- plot.pcs.from.betas(pcs, samplesheet,variables, cols)
    colnames(res)[which(colnames(res) == "x")] <- "batch.variable"
    colnames(res)[which(colnames(res) == "l")] <- "batch.value"
    colnames(res)[which(colnames(res) == "y")] <- "principal.component"
    ret$tab <- res
    ret
}


plot.pcs.from.betas <- function(pcs, samplesheet,variables,cols=NULL) {
    stopifnot(is.matrix(pcs) && nrow(pcs) == nrow(samplesheet))
    stopifnot(is.numeric(pcs))
    stopifnot(samplesheet$Sample_Name==names(pcs))
   
    dat<-samplesheet[,which(names(samplesheet)%in%variables)]
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

#' Specify parameters for testing normalization
#'
#' @param  batch.threshold Default = 1e-50. Which pvalue threshold to show in table
#' @param  batch.pcs Default = 1:10. Number of PCs to test against batch variables
#' @param  colours Colours to use for scatterplots.
#' @export
#' @return List of parameters
#' @examples \dontrun{
#'
#'}
meffil.normalization.parameters.from.betas <- function(
                                            batch.pcs = 1:10,
                                            batch.threshold=1e-50,
                                            colours=NULL) {                                           
    
    parameters <- list(
        batch.pcs = batch.pcs,
        batch.threshold=batch.threshold,
        colours=colours
  )
    return(parameters)
}
