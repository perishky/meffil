#' @export
meffil.ewas.report <- function(ewas.summary,
                               output.file = "ewas-report.html",
                               author = "Analyst",
                               study = "IlluminaHuman450 data",
                               ...) {
    meffil:::msg("Writing report as html file to", output.file)
    path <- system.file("reports", package="meffil")
    knit.report(file.path(path, "ewas-report.rmd"),output.file, ...)
}

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


#' @export
meffil.ewas.parameters <- function(sig.threshold=1e-7,max.plots=10) {
    list(sig.threshold=sig.threshold,
         max.plots=max.plots)
}
