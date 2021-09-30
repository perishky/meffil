#' Epigenome-wide association study
#'
#' Test association with each CpG site.
#'
#' @param beta Methylation levels matrix,
#' one row per CpG site, one column per sample
#' or the filename of GDS (Genomic Data Structure) output from
#' \code{\link{meffil.normalize.samples}}.
#' @param variable Independent variable vector.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param batch Batch vector to be included as a random effect (Default: NULL). Ignored if \code{beta} is a GDS filename.
#' @param weights Non-negative observation weights.
#' Can be a numeric matrix of individual weights of same dimension as \code{beta},
#' or a numeric vector of weights with length \code{ncol(beta)},
#' or a numeric vector of weights with length \code{nrow(beta)}. 
#' @param sites Restrict the EWAS to the given CpG sites -- must match row names of \code{beta} (Default: NULL).
#' @param samples Restrict the EWAS to the given samples -- must match column names of \code{beta} (Default: NULL). 
#' @param cell.counts Proportion of cell counts for one cell type in cases
#' where the samples are mainly composed of two cell types (e.g. saliva) (Default: NULL). Ignored if \code{beta} is a GDS filename.
#' @param isva Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: FALSE).  
#' @param sva Apply Surrogate Variable Analysis (SVA) to the
#' methylation levels and covariates and include
#' the resulting variables as covariates in a regression model (Default: TRUE).
#' @param smartsva Apply the SmartSVA algorithm to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: FALSE).
#' @param smartsva.alpha alpha argument to SmartSVA providing the
#' initial point for optimization.  Smaller values
#' reduce the number of iterations needed to reach convergence.
#' Setting this 1 will produce exactly the outputs as SVA. (Default: 0.5).
#' @param n.sv Number of surrogate variables to calculate (Default: NULL).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05).
#' @param robust Test associations with the 'robust' option when \code{\link{limma::eBayes}}
#' is called (Default: TRUE). Ignored if \code{beta} is a GDS filename.
#' @param rlm If \code{beta} is a matrix, then test associations with
#' the 'robust' option when \code{\link{limma:lmFit}} is called.
#' If \code{beta} is a GDS filename, then test associations
#' using robust regression using \code{\link{MASS::rlm}}
#' and calculate statistical significance using \code{\link{lmtest::coeftest}}
#' with \code{vcov=\link{sandwich::vcovHC}(fit, type="HC0")}
#' (Default: FALSE).
#' @param outlier.iqr.factor For each CpG site, prior to fitting regression models,
#' set methylation levels less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA).
#' @param most.variable Apply (Independent) Surrogate Variable Analysis to the 
#' given most variable CpG sites (Default: 50000).
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}  (Default: NA).
#' @param verbose Set to TRUE if status updates to be printed (Default: FALSE).
#'
#' @export
meffil.ewas <- function(beta, variable,
                        covariates=NULL, batch=NULL, weights=NULL,
                        sites=NULL, samples=NULL,
                        cell.counts=NULL,
                        isva=F, sva=T, smartsva=F, ## cate?
                        smartsva.alpha=0.5,
                        n.sv=NULL,
                        winsorize.pct=0.05,
                        robust=FALSE,
                        rlm=FALSE,
                        outlier.iqr.factor=NA, ## typical value = 3
                        most.variable=50000,
                        featureset=NA,
                        random.seed=20161123,
                        lmfit.safer=F,
                        verbose=F) {

    if (is.matrix(beta)) {
        all.sites <- rownames(beta)
        all.samples <- colnames(beta)
    }
    else {
        beta.dims <- meffil:::retrieve.gds.dims(beta)
        all.sites <- beta.dims[[1]]
        all.samples <- beta.dims[[2]]

        if (!is.null(batch))
            stop("'batch' is ignored when 'beta' is a GDS filename")
        if (!is.null(cell.counts))
            stop("'cell.counts' is ignored when 'beta' is a GDS filename")
        if (robust)
            warning("'robust' is ignored when 'beta' is a GDS filename")
    }
    
    if (is.na(featureset))
        featureset <- guess.featureset(all.sites)
    features <- meffil.get.features(featureset)

    if (!is.null(sites))
        sites <- intersect(all.sites, sites)
    else
        sites <- all.sites

    stopifnot(length(sites) > 1 && all(sites %in% features$name))
    
    if (is.null(samples))
        samples <- all.samples
        
    stopifnot(length(samples) > 1 && all(samples %in% all.samples))

    stopifnot(length(all.samples) == length(variable))
    stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == length(all.samples))
    stopifnot(is.null(weights)
              || (is.numeric(weights)
                  && is.vector(weights)
                  && length(weights) == length(all.samples)))
    stopifnot(is.null(batch) || length(batch) == length(all.samples))
    stopifnot(most.variable > 1)
    stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && winsorize.pct < 0.5)
    
    original.variable <- variable
    original.covariates <- covariates
    if (is.character(variable))
        variable <- as.factor(variable)
    
    stopifnot(!is.factor(variable) || is.ordered(variable) || length(levels(variable)) == 2)
    
    meffil:::msg("Simplifying any categorical variables.", verbose=verbose)
    variable <- meffil:::simplify.variable(variable)
    if (!is.null(covariates))
        covariates <- do.call(cbind, lapply(covariates, meffil:::simplify.variable))

    sample.idx <- which(!is.na(variable) & all.samples %in% samples)
    if (!is.null(covariates))
        sample.idx <- intersect(sample.idx, which(apply(!is.na(covariates), 1, all)))
    
    meffil:::msg("Removing", length(samples) - length(sample.idx), "missing case(s).", verbose=verbose)
    weights <- weights[sample.idx]    
    samples <- all.samples[sample.idx]
    variable <- variable[sample.idx]
    original.variable <- original.variable[sample.idx]
    original.covariates <- original.covariates[sample.idx,,drop=F]

    if (!is.null(covariates))
        covariates <- covariates[sample.idx,,drop=F]

    if (!is.null(batch))
        batch <- batch[sample.idx]

    if (!is.null(cell.counts))
        cell.counts <- cell.counts[sample.idx]

    if (!is.null(covariates)) {
        pos.var.idx <- which(apply(covariates, 2, var, na.rm=T) > 0)
        meffil:::msg("Removing", ncol(covariates) - length(pos.var.idx), "covariates with no variance.",
            verbose=verbose)
        covariates <- covariates[,pos.var.idx, drop=F]
    }

    covariate.sets <- list(none=NULL)
    if (!is.null(covariates))
        covariate.sets$all <- covariates
    
    surrogates.ret <- NULL
    if (isva || sva || smartsva) {
        meffil:::msg("Computing surrogate variables ...", verbose=verbose)
        sva.sites <- intersect(sites, meffil:::autosomal.sites(beta))
        sva.sets <- add.surrogate.variables(
            beta, sva.sites, samples,
            most.variable,
            variable, covariates,
            winsorize.pct, outlier.iqr.factor,
            random.seed, isva, sva, smartsva, n.sv,
            smartsva.alpha, verbose)
        surrogates.ret <- sva.sets$last
        covariate.sets <- c(covariate.sets, sva.sets$sets)
    }

    msg("Starting EWAS ...", verbose=verbose)

    if (is.matrix(beta)) {
        analyses <- meffil:::ewas.by.matrix(
            variable, beta, sites, samples,
            covariate.sets,
            batch, weights, cell.counts,
            winsorize.pct, outlier.iqr.factor,
            robust, rlm, lmfit.safer, verbose)
    }
    else {
        analyses <- meffil:::ewas.by.gds(
            variable, beta, sites, samples,
            covariate.sets, weights,
            winsorize.pct, outlier.iqr.factor,
            rlm, verbose)
    }

    msg("Finished EWAS.", verbose=verbose)

    p.values <- sapply(analyses, function(analysis) analysis$table$p.value)
    coefficients <- sapply(analyses, function(analysis) analysis$table$coefficient)
    rownames(p.values) <- rownames(coefficients) <- rownames(analyses[[1]]$table)

    for (name in names(analyses)) {
        idx <- match(rownames(analyses[[name]]$table), features$name)
        analyses[[name]]$table$chromosome <- features$chromosome[idx]
        analyses[[name]]$table$position <- features$position[idx]
    }

    list(class="ewas",
         version=packageVersion("meffil"),
         samples=samples,
         variable=original.variable,
         covariates=original.covariates,
         winsorize.pct=winsorize.pct,
         batch=batch,
         robust=robust,
         rlm=rlm,
         weights=weights,
         outlier.iqr.factor=outlier.iqr.factor,
         most.variable=most.variable,
         p.value=p.values,
         coefficient=coefficients,
         analyses=analyses,
         random.seed=random.seed,
         sva.ret=surrogates.ret)
}

is.ewas.object <- function(object)
    is.list(object) && all(c("class","version","samples","variable","covariates",
                             "winsorize.pct","outlier.iqr.factor",
                             "analyses","random.seed")
                           %in% names(object))

ewas.by.matrix <- function(variable, beta, sites, samples,
                           covariate.sets, batch, weights, cell.counts,
                           winsorize.pct, outlier.iqr.factor,
                           robust, rlm, lmfit.safer, verbose) {

    stopifnot(all(!is.na(variable)))
    stopifnot(length(variable) == length(samples))
    stopifnot(is.null(batch) || length(batch) == length(samples))
    stopifnot(is.null(weights) || length(weights) == length(variable))
    stopifnot(is.null(cell.counts)
              || length(cell.counts) == length(samples)
              && all(cell.counts >= 0 & cell.counts <= 1))
    stopifnot(length(sites) > 1 && all(sites %in% rownames(beta)))
    stopifnot(length(samples) > 1 && all(samples %in% colnames(beta)))
    stopifnot(all(sapply(covariate.sets, function(set) is.null(set) || nrow(set) == length(samples))))
    stopifnot(is.logical(robust))
    stopifnot(is.logical(rlm))
    stopifnot(is.logical(lmfit.safer))

    beta <- beta[sites,samples]
    
    if (is.numeric(winsorize.pct) || is.numeric(outlier.iqr.factor)) {
        meffil:::msg("Handling methylation outliers", verbose=verbose) 
        beta <- meffil.handle.outliers(
            beta,
            winsorize.pct,
            outlier.iqr.factor)       
    }
    
    sapply(names(covariate.sets), function(name) {
        meffil:::msg("EWAS for covariate set", name, verbose=verbose)
        covariates <- covariate.sets[[name]]
        meffil:::ewas.limma(
            variable,
            beta=beta,
            covariates=covariates,
            batch=batch,
            weights=weights,
            cell.counts=cell.counts,
            winsorize.pct=winsorize.pct, 
            robust=robust, 
            rlm=rlm, 
            lmfit.safer=lmfit.safer,
            verbose=verbose)
    }, simplify=F)
}


# Test associations between \code{variable} and each row of \code{beta}
# while adjusting for \code{covariates} (fixed effects) and \code{batch} (random effect).
# If \code{cell.counts} is not \code{NULL}, then it is assumed that
# the methylation data is derived from samples with two cell types.
# \code{cell.counts} should then be a vector of numbers
# between 0 and 1 of length equal to \code{variable} corresponding
# to the proportions of cell of a selected cell type in each sample.
# The regression model is then modified in order to identify
# associations specifically in the selected cell type (PMID: 24000956).
ewas.limma <- function(variable, beta, covariates,
                       batch, weights, cell.counts,
                       winsorize.pct,
                       robust, rlm, lmfit.safer, verbose) {
    stopifnot(all(!is.na(variable)))
    stopifnot(length(variable) == ncol(beta))
    stopifnot(is.null(covariates) || nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(is.null(cell.counts)
              || length(cell.counts) == ncol(beta)
              && all(cell.counts >= 0 & cell.counts <= 1))
  
    method <- "ls"
    if (rlm) method <- "robust"
  
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
        
        design <- cbind(design * cell.counts,
                        typeB=design * (1-cell.counts))
    }

    meffil:::msg("Linear regression with limma::lmFit", verbose=verbose)
    fit <- NULL
    batch.cor <- NULL
    if (!is.null(batch)) {
        meffil:::msg("Adjusting for batch effect", verbose=verbose)
        corfit <- duplicateCorrelation(beta, design, block=batch, ndups=1)
        batch.cor <- corfit$consensus
        meffil:::msg("Linear regression with batch as random effect", verbose=verbose)
        tryCatch({
            fit <- lmFit(
                beta, design, method=method,
                block=batch, cor=batch.cor, weights=weights)
        }, error=function(e) {
            print(e)
            meffil:::msg("lmFit failed with random effect batch variable, omitting", verbose=verbose)
        })
    }
    if (is.null(fit)) {
        meffil:::msg("Linear regression with only fixed effects", verbose=verbose)
        batch <- NULL
        if (!lmfit.safer)
           fit <- lmFit(beta, design, method=method, weights=weights)
        else
           fit <- lmfit.safer(beta, design, method=method, weights=weights, verbose=verbose)
    }
    
    meffil:::msg("Empirical Bayes", verbose=verbose)
    if (is.numeric(winsorize.pct) && robust) {
      fit.ebayes <- eBayes(fit, robust=T, winsor.tail.p=c(winsorize.pct, winsorize.pct))
    } else {
      fit.ebayes <- eBayes(fit, robust=robust)
    }

    alpha <- 0.975
    std.error <- (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,"variable"])
    margin.error <- (std.error * qt(alpha, df=fit.ebayes$df.total))
    n <- apply(beta, 1, function(v) sum(!is.na(v)))

    list(design=design,
         batch=batch,
         batch.cor=batch.cor,
         cell.counts=cell.counts,
         table=data.frame(p.value=fit.ebayes$p.value[,"variable"],
             fdr=p.adjust(fit.ebayes$p.value[,"variable"], "fdr"),
             p.holm=p.adjust(fit.ebayes$p.value[,"variable"], "holm"),
             t.statistic=fit.ebayes$t[,"variable"],
             coefficient=fit.ebayes$coefficient[,"variable"],
             coefficient.ci.high=fit.ebayes$coefficient[,"variable"] + margin.error,
             coefficient.ci.low=fit.ebayes$coefficient[,"variable"] - margin.error,
             coefficient.se=std.error,
             n=n))
}

## apply lmFit to subsets of the methylation matrix to hopefully avoid 
## out-of-memory errors (mainly for really large datasets or low memory servers)
lmfit.safer <- function(beta, design, method, weights, n.partitions=8, verbose=F) {
  partitions <- sample(1:n.partitions, nrow(beta), replace=T)
  fits <- lapply(1:n.partitions, function(part) {
    meffil:::msg("Applying limma::lmFit to partition ", part, " of ", n.partitions, ".", verbose=verbose)
    idx <- which(partitions == part)
    lmFit(beta[idx,], design, method=method, weights=weights)
  })
  idx <- unlist(lapply(1:n.partitions, function(part) which(partitions==part)))
  idx <- match(1:nrow(beta), idx)
  fit <- fits[[1]]
  for (item in c("coefficients", "stdev.unscaled")) {
      fit[[item]] <- do.call(rbind, lapply(fits, function(fit) fit[[item]]))[idx,]
  }
  for (item in c("df.residual", "sigma", "Amean")) {
     fit[[item]] <- do.call(c, lapply(fits, function(fit) fit[[item]]))[idx]
  }
  fit
}


# Test associations between \code{variable} and each row of \code{beta}
# while adjusting for \code{covariates}.
ewas.by.gds <- function(variable, beta, sites, samples,
                        covariate.sets, weights,
                        winsorize.pct, outlier.iqr.factor,
                        rlm, verbose) { 
    stopifnot(all(!is.na(variable)))

    stopifnot(length(variable) == length(samples))
    stopifnot(all(sapply(covariate.sets, function(set) is.null(set) || nrow(set) == length(samples))))
    stopifnot(is.null(weights) || length(weights) == length(variable))
    
    designs <- lapply(covariate.sets, function(covariates) {
        if (is.null(covariates))
            design <- data.frame(intercept=1, variable=variable)
        else
            design <- data.frame(intercept=1, variable=variable, covariates)
        rownames(design) <- samples
        design
    })        

    stats <- meffil:::lapply.gds(
        beta, margin=1, sites=sites, samples=samples,
        type="list",
        FUN=ewas.test.models,        
        AFUN=ifelse(rlm,meffil:::ewas.test.rlm,meffil:::ewas.test.glm),
        designs=designs,
        weights=weights,
        winsorize.pct=winsorize.pct,
        outlier.iqr.factor=outlier.iqr.factor)
    stats <- do.call(rbind, stats)
    
    analyses <- lapply(1:length(designs), function(i) {
        idx <- seq(i,nrow(stats),length(designs))
        stats <- data.frame(stats[idx,])
        rownames(stats) <- sites
        stats$fdr <- p.adjust(stats$p.value, "fdr")
        stats$p.holm <- p.adjust(stats$p.value, "holm")                
        list(design=designs[[i]],
             table=stats)
    })
    names(analyses) <- names(covariate.sets)
    analyses
}

ewas.test.models <- function(cpg,AFUN,designs,weights,winsorize.pct,outlier.iqr.factor,...) {
    cpg <- meffil::meffil.handle.outliers(cpg,winsorize.pct, outlier.iqr.factor)
    idx <- which(!is.na(cpg))
    cpg <- cpg[idx]
    if (!is.null(weights)) weights <- weights[idx]
    t(sapply(designs, function(design) {
        AFUN(cpg,design[idx,],weights,...)
    }))
}

ewas.test.glm <- function(cpg,design,weights,...) {
    fit <- glm.fit(x=design,y=cpg,weights=weights,...)
    ret <- coef(summary.glm(fit))["variable",]
    names(ret) <- c("coefficient","coefficient.se","t.statistic","p.value")
    c(ret,
      coefficient.ci.low=ret[["coefficient"]]-1.96*ret[["coefficient.se"]],
      coefficient.ci.high=ret[["coefficient"]]+1.96*ret[["coefficient.se"]],
      n=sum(!is.na(cpg)))
}

ewas.test.rlm <- function(cpg,design,weights=NULL,...) {
    tryCatch({
        fit <- MASS::rlm(x=design,y=cpg,...)
        ret <- lmtest::coeftest(fit, vcov=sandwich::vcovHC(fit, type="HC0"))
        ci <- confint(ret)["variable",]
        ret <- ret["variable",]
        names(ret) <- c("coefficient","coefficient.se","t.statistic","p.value") 
        c(ret,
          coefficient.ci.low=ci[["2.5 %"]],
          coefficient.ci.high=ci[["97.5 %"]],
          n=sum(!is.na(cpg)))
    }, error=function(e) {
        c("coefficient"=NA,"coefficient.se"=NA,"t.statistic"=NA,"p.value"=NA,
          "coefficient.ci.low"=NA,"coefficient.ci.high"=NA,n=sum(!is.na(cpg)))
    })
}
