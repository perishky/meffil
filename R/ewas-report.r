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
    path <- system.file("reports", package="meffil")
    knit.report(file.path(path, "ewas-report.rmd"),output.file, ...)
}

#' Summarize EWAS results.
#'
#' Generates variable and covariate summary tables, QQ plots, Manhattan plots, a list of associations,
#' plots of the strongest associations and plots of selected CpG sites.
#'
#' @param ewas.object From \code{\link{meffil.ewas}()}.
#' @param selected.cpg.sites Vector of CpG site names to plot (Default: character(0)).
#' @param  parameters Default = meffil.ewas.parameters(). List of parameter values. See \code{\link{meffil.ewas.parameters}()}.
#' @export
#' @return List
#' 
#' @export
meffil.ewas.summary <- function(ewas.object, beta,
                                selected.cpg.sites=character(0),
                                parameters=meffil.ewas.parameters(),
                                verbose=T) {
    min.p <- rowMins(ewas.object$p.value)
    parameters$practical.threshold <- min.p[order(min.p, decreasing=F)][parameters$max.plots+1]
    if (parameters$practical.threshold > parameters$sig.threshold) {
        parameters$practical.threshold <- sig.threshold
    }
        
    sig.idx <- unique(which(ewas.object$p.value < parameters$sig.threshold, arr.ind=T)[,"row"])
    sig.cpg.stats <- list(p.value=ewas.object$p.value[sig.idx,],
                          coefficient=ewas.object$coefficient[sig.idx,])

    selected.cpg.idx <- match(selected.cpg.sites, rownames(ewas.object$p.value))
    selected.cpg.stats <- list(p.value=ewas.object$p.value[selected.cpg.idx,],
                               coefficient=ewas.object$coefficient[selected.cpg.idx,])

    msg("QQ plots", verbose=T)
    qq.plots <- meffil.ewas.qq.plot(ewas.object,
                                    sig.threshold=parameters$sig.threshold)

    msg("Manhattan plots", verbose=T)
    manhattan.plots <- meffil.ewas.manhattan.plot(ewas.object,
                                                  sig.threshold=parameters$sig.threshold)

    practical.idx <- which(rowMins(sig.cpg.stats$p.value) < parameters$practical.threshold)
    cpg.sites <- union(rownames(sig.cpg.stats$p.value)[practical.idx], 
                       rownames(selected.cpg.stats$p.value))
    msg("CpG site plots:", length(cpg.sites), verbose=T)
    cpg.plots <- sapply(cpg.sites, function(cpg) {
        msg("Plotting", cpg, verbose=T)
        meffil.ewas.cpg.plot(ewas.object, cpg, beta=beta)
    }, simplify=F)

    msg("Sample characteristics", verbose=T)
    sample.characteristics <- meffil.ewas.sample.characteristics(ewas.object)
    covariate.associations <- meffil.ewas.covariate.associations(ewas.object)

    cpg.sites <- union(rownames(sig.cpg.stats$p.value),
                       rownames(selected.cpg.stats$p.value))
    cpg.sites <- ewas.object$analyses[[1]]$table[cpg.sites,c("chromosome","position")]

    parameters$winsorize.pct <- ewas.object$winsorize.pct
    
    list(parameters=parameters,
         qq.plots=qq.plots,
         manhattan.plots=manhattan.plots,
         sig.cpg.stats=sig.cpg.stats,
         selected.cpg.stats=selected.cpg.stats,
         cpg.plots=cpg.plots,
         cpg.sites=cpg.sites,
         sample.characteristics=sample.characteristics,
         covariate.associations=covariate.associations)         
}

#' Specify parameters for QC
#'
#' 
#' @param sig.threshold P-value threshold for significance (Default: 1e-7).
#' @param max.plots Maximum number of plots to generate (Default: 10).
#' @return List of parameter values
#'
#' @export
meffil.ewas.parameters <- function(sig.threshold=1e-7,max.plots=10) {
    list(sig.threshold=sig.threshold,
         max.plots=max.plots)
}
