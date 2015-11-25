#' Epigenome-wide association study
#'
#' Test association with each CpG site.
#'
#' @param beta Methylation levels matrix, one row per CpG site, one column per sample.
#' @param variable Independent variable vector.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param batch Batch vector to be included as a random effect (Default: NULL).
#' @param save.beta Number of most strongly associated CpG sites for which to retain
#' methylation levels for later plotting (Default: 100).
#' @param cell.counts Proportion of cell counts for one cell type in cases
#' where the samples are mainly composed of two cell types (e.g. saliva) (Default: NULL).
#' @param isva0 Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: TRUE).
#' @param isva1 Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels, covariates and \code{isva0} variables and include
#' the resulting variables as covariates in a regression model (Default: TRUE).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level (Default: 0.05).
#' @param most.variable Apply Independent Surrogate Variable Analysis to the 
#' given most variable CpG sites (Default: 50000).
#'
#' @export
meffil.ewas <- function(beta, variable,
                        covariates=NULL, batch=NULL,
                        save.beta=100, cell.counts=NULL,
                        isva0=T, isva1=T,
                        winsorize.pct=0.05, ## perhaps better, winsorize at 25-percentile - iqr?
                        most.variable=min(nrow(beta), 50000),
                        featureset=NULL,
                        verbose=F) {

    if (is.null(featureset))
        featureset <- guess.architecture(beta)
    features <- meffil.get.features(featureset)
    
    stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% features$name))
    stopifnot(ncol(beta) == length(variable))
    stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(most.variable > 1 && most.variable <= nrow(beta))
    stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && winsorize.pct < 0.5)

    original.variable <- variable
    original.covariates <- covariates
    if (is.character(variable))
        variable <- as.factor(variable)
    
    stopifnot(!is.factor(variable) || is.ordered(variable) || length(levels(variable)) == 2)
    
    simplify.variable <- function(v) {
        if (is.character(v))
            v <- as.factor(v)
        if (is.factor(v)) {
            v <- droplevels(v)
            if (is.ordered(v) || length(levels(v)) <= 2)
                as.integer(v) - 1
            else
                model.matrix(~ 0 + v)[,-1,drop=F]
        } else
            v
    }

    msg("Simplifying any categorical variables.", verbose=verbose)
    variable <- simplify.variable(variable)
    if (!is.null(covariates))
        covariates <- do.call(cbind, lapply(covariates, simplify.variable))

    sample.idx <- which(!is.na(variable))
    msg("Removing", ncol(beta) - length(sample.idx), "missing case(s).", verbose=verbose)
    beta <- beta[,sample.idx]
    variable <- variable[sample.idx]

    if (!is.null(covariates))
        covariates <- covariates[sample.idx,]

    if (!is.null(batch))
        batch <- batch[sample.idx]

    if (!is.null(cell.counts))
        cell.counts <- cell.counts[sample.idx]

    if (!is.null(covariates)) {
        pos.var.idx <- which(apply(covariates, 2, var) > 0)
        msg("Removing", ncol(covariates) - length(pos.var.idx), "covariates with no variance.",
            verbose=verbose)
        covariates <- covariates[,pos.var.idx, drop=F]
    }

    covariate.sets <- list(none=NULL)

    if (is.numeric(winsorize.pct))  {
        msg(winsorize.pct, "- winsorizing the beta matrix.", verbose=verbose)
        beta <- winsorize(beta, pct=winsorize.pct)
    }

    if (isva0 || isva1) {
        msg("ISVA with no covariates.", verbose=verbose)
        var.idx <- order(rowVars(beta, na.rm=T), decreasing=T)[1:most.variable]    
        isva0 <- DoISVA(beta[var.idx,,drop=F], variable, verbose=verbose)
        covariate.sets$isva0 <- as.data.frame(isva0$isv)
    }

    if (isva1) {
        msg("ISVA with covariates.", verbose=verbose)
        if (!is.null(covariates)) {
            factor.log <- rep(FALSE, nrow(covariates))
            isva1 <- DoISVA(beta[var.idx,,drop=F], variable, cf.m=covariates, factor.log=factor.log,
                        verbose=verbose)
            covariate.sets$isva1 <- as.data.frame(isva1$isv)
        }
    }

    if (!is.null(covariates))
        covariate.sets$all <- covariates

    analyses <- sapply(names(covariate.sets), function(name) {
        msg("EWAS for covariate set", name, verbose=verbose)
        covariates <- covariate.sets[[name]]
        ewas(variable,
             beta=beta,
             covariates=covariates,
             batch=batch,
             save.beta=save.beta,
             cell.counts=cell.counts)
    }, simplify=F)

    p.values <- sapply(analyses, function(analysis) analysis$table$p.value)
    coefficients <- sapply(analyses, function(analysis) analysis$table$coefficient)
    rownames(p.values) <- rownames(coefficients) <- rownames(analyses[[1]]$table)

    for (name in names(analyses)) {
        idx <- match(rownames(analyses[[name]]$table), features$name)
        analyses[[name]]$table$chromosome <- features$chromosome[idx]
        analyses[[name]]$table$position <- features$position[idx]
    }
    
    list(class="ewas",
         samples=sample.idx,
         variable=original.variable[sample.idx],
         covariates=original.covariates[sample.idx,],
         winsorize.pct=winsorize.pct,
         most.variable=most.variable,
         p.value=p.values,
         coefficient=coefficients,
         analyses=analyses)
}

is.ewas.object <- function(object)
    is.list(object) && "class" %in% names(object) && object$class == "ewas"


ewas <- function(variable, beta, covariates=NULL, batch=NULL, save.beta=100, cell.counts=NULL,
                 verbose=F) {
    stopifnot(all(!is.na(variable)))
    stopifnot(length(variable) == ncol(beta))
    stopifnot(is.null(covariates) || nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(is.null(cell.counts)
              || length(cell.counts) == ncol(beta)
              && all(cell.counts >= 0 & cell.counts <= 1))
    stopifnot(save.beta > 0)

    if (is.null(covariates))
        design <- data.frame(intercept=1, variable=variable)
    else
        design <- data.frame(intercept=1, variable=variable, covariates)
    rownames(design) <- colnames(beta)

    if (!is.null(cell.counts)) {
        ## Irizarray method: Measuring cell-type specific differential methylation
        ##      in human brain tissue
        ## Mi is methylation level for sample i
        ## Xi is the value of the variable of interest for sample i
        ## pi is the proportion of the target cell type for sample i
        ## thus: Mi ~ target cell methylation * pi + other cell methylation * (1-pi)
        ## and target cell methylation = base target methylation + effect of Xi value
        ## and other cell methylation = base other methylation + effect of Xi value
        ## thus:
        ## Mi = (T + U Xi)pi + (O + P Xi)(1-pi) + e
        ##    = C + A Xi pi + B Xi (1-pi) + e        
        design <- design[,-which(colnames(design) == "intercept")]
        
        design <- cbind(typeA=design * cell.counts,
                        typeB=design * (1-cell.counts))
    }

    msg("Linear regression with limma::lmFit", verbose=verbose)
    fit <- NULL
    batch.cor <- NULL
    if (!is.null(batch)) {
        msg("Adjusting for batch effect", verbose=verbose)
        corfit <- duplicateCorrelation(beta, design, block=batch, ndups=1)
        batch.cor <- corfit$consensus
        msg("Linear regression with batch as random effect", verbose=verbose)
        tryCatch(fit <- lmFit(beta, design, block=batch, cor=batch.cor),
                 error=function(e) {
                     print(e)
                     msg("lmFit failed with random effect batch variable, omitting", verbose=verbose)
                 })
    }
    if (is.null(fit)) {
        msg("Linear regression with only fixed effects", verbose=verbose)
        batch <- NULL
        fit <- lmFit(beta, design)
    }
    
    msg("Empirical Bayes", verbose=verbose)
    fit.ebayes <- eBayes(fit, robust=T)

    alpha <- 0.975
    margin.error <- (sqrt(fit.ebayes$s2.post)
                     * fit.ebayes$stdev.unscaled[,"variable"]
                     * qt(alpha, df=fit.ebayes$df.total))

    save.idx <- order(fit.ebayes$p.value[,"variable"], decreasing=F)[1:save.beta]    
    
    list(design=design,
         batch=batch,
         batch.cor=batch.cor,
         cell.counts=cell.counts,
         beta=beta[save.idx,,drop=F],
         table=data.frame(p.value=fit.ebayes$p.value[,"variable"],
             fdr=p.adjust(fit.ebayes$p.value[,"variable"], "fdr"),
             p.holm=p.adjust(fit.ebayes$p.value[,"variable"], "holm"),
             t.statistic=fit.ebayes$t[,"variable"],
             coefficient=fit.ebayes$coefficient[,"variable"],
             coefficient.ci.high=fit.ebayes$coefficient[,"variable"] + margin.error,
             coefficient.ci.low=fit.ebayes$coefficient[,"variable"] - margin.error))
}




