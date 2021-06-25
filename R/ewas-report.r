#' Generate EWAS report.
#'
#' Generate HTML file that summarises the EWAS. 
#'
#' @param  ewas.summary Output from \code{meffil.ewas.summary}.
#' @param  output.file Default = "ewas-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "Illumina methylation data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' 
#' @export
meffil.ewas.report <- function(ewas.summary,
                               output.file = "ewas-report.html",
                               author = "Analyst",
                               study = "Illumina methylation data",
                               ...) {
    meffil:::msg("Writing report as html file to", output.file)
    report.path <- system.file("reports", package="meffil")
    require(knitr)
    require(Cairo)
    require(gridExtra)
    opts <- opts_chunk$get()
    on.exit(opts_chunk$set(opts))
    opts_chunk$set(warning=FALSE, echo=FALSE, message=FALSE, results="asis", fig.width=6, fig.height=6, dev="CairoPNG")
    knit.report(file.path(report.path, "ewas-report.rmd"),output.file, ...)
}

#' Summarize EWAS results.
#'
#' Generates variable and covariate summary tables, QQ plots, Manhattan plots, a list of associations,
#' plots of the strongest associations and plots of selected CpG sites.
#'
#' @param ewas.object From \code{\link{meffil.ewas}()}.
#' @param beta Methylation levels used in the analysis,
#' either a matrix with
#' one row per CpG site and one column per sample
#' or the filename of a GDS file (Genomic Data Structure).
#' @param selected.cpg.sites Vector of CpG site names to plot (Default: character(0)).
#' @param parameters Default = meffil.ewas.parameters(). List of parameter values. See \code{\link{meffil.ewas.parameters}()}.
#' @export
#' @return List
#' 
#' @export
meffil.ewas.summary <- function(ewas.object, beta,
                                selected.cpg.sites=character(0),
                                parameters=meffil.ewas.parameters(),
                                verbose=T) {
    stopifnot(parameters$model %in% colnames(ewas.object$p.value))
    stopifnot(parameters$max.plots < nrow(ewas.object$p.value))
    stopifnot(all(selected.cpg.sites %in% rownames(ewas.object$p.value)))
    
    p.values <- ewas.object$p.value[,parameters$model]
    parameters$practical.threshold <- p.values[order(p.values)[parameters$max.plots+1]]

    if (is.na(parameters$sig.threshold)) 
        parameters$sig.threshold <- 0.05/nrow(ewas.object$p.value)
        
    sig.idx <- unique(which(ewas.object$p.value < parameters$sig.threshold, arr.ind=T)[,"row"])
    practical.idx <- which(p.values < parameters$practical.threshold)
    selected.idx <- match(selected.cpg.sites, rownames(ewas.object$p.value))    

    significant.sites <- rownames(ewas.object$p.value)[sig.idx]
    selected.sites <- rownames(ewas.object$p.value)[selected.idx]
  
    cpg.idx <- union(sig.idx, union(practical.idx, selected.idx))
    cpg.sites <- rownames(ewas.object$p.value)[cpg.idx]

    cpg.stats <- data.frame(ewas.object$analyses[[1]]$table[cpg.sites, c("chromosome","position")],
                            p.value=ewas.object$p.value[cpg.idx,],
                            coefficient=ewas.object$coefficient[cpg.idx,])
    
    msg("QQ plots", verbose=T)
    qq.plots <- meffil.ewas.qq.plot(ewas.object,
                                    sig.threshold=parameters$sig.threshold,
                                    lambda.method=parameters$qq.inflation.method)

    msg("Manhattan plots", verbose=T)
    manhattan.plots <- meffil.ewas.manhattan.plot(ewas.object,
                                                  sig.threshold=parameters$sig.threshold)

    plot.sites <- rownames(ewas.object$p.value)[union(practical.idx, selected.idx)]
    msg("CpG site plots:", length(plot.sites), verbose=T)
    if (is.character(beta)) 
        beta <- retrieve.gds.methylation(beta, sites=plot.sites, samples=NULL)
    cpg.plots <- sapply(plot.sites, function(cpg) {
        msg("Plotting", cpg, verbose=T)
        meffil.ewas.cpg.plot(ewas.object, cpg=cpg, beta=beta)
    }, simplify=F)

    msg("Sample characteristics", verbose=T)
    sample.characteristics <- meffil.ewas.sample.characteristics(ewas.object)
    covariate.associations <- meffil.ewas.covariate.associations(ewas.object)

    parameters$winsorize.pct <- ewas.object$winsorize.pct
    parameters$outlier.iqr.factor <- ewas.object$outlier.iqr.factor
    parameters$rlm <- ewas.object$rlm
    parameters$most.variable <- ewas.object$most.variable
    parameters$random.seed <- ewas.object$random.seed
    parameters$sample.size <- length(ewas.object$samples)
   
    ## order CpG sites and plots by the default model p-value
    sort.by.p <- function(x) {
      sites <- x
      if (!is.character(x))
        sites <- names(x)
      idx <- match(sites, rownames(ewas.object$p.value))
      p <- ewas.object$p.value[idx,parameters$model]
      x[order(p)]
    }
  
    list(parameters=parameters,
         qq.plots=qq.plots,
         manhattan.plots=manhattan.plots,
         cpg.stats=cpg.stats,
         cpg.plots=sort.by.p(cpg.plots),
         significant.sites=sort.by.p(significant.sites),
         selected.sites=sort.by.p(selected.sites),         
         sample.characteristics=sample.characteristics,
         covariate.associations=covariate.associations)         
}

#' Specify parameters for QC
#'
#' 
#' @param sig.threshold P-value threshold for significance (Default: NA).
#' If NA, then threshold used will be 0.05 divided by the number of tests/probes.
#' @param max.plots Maximum number of plots to generate (Default: 10).
#' @param model Model to use for selecting associations: "none" (no covariates),
#' "all" (all covariates), "isva" (independent surrogate variables),
#' and "sva" (surrogate variables) (Default: "none").
#' @param qq.inflation.method Method for calculating genomic inflation lambda.
#' Valid values are "median", "regression" or "robust" (Default: "median").
#' @return List of parameter values
#'
#' @export
meffil.ewas.parameters <- function(sig.threshold=NA,
                                   max.plots=10,
                                   model="none",
                                   qq.inflation.method="median") {
    list(sig.threshold=sig.threshold,
         max.plots=max.plots,
         model=model,
         qq.inflation.method=qq.inflation.method)
}
