#' Epigenome-wide association study
#'
#' Test association with each CpG site.
#'
#' @param beta Methylation levels matrix, one row per CpG site, one column per sample.
#' @param variable Independent variable vector.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param batch Batch vector to be included as a random effect (Default: NULL).
#' @param weights Non-negative observation weights.
#' Can be a numeric matrix of individual weights of same dimension as \code{beta},
#' or a numeric vector of weights with length \code{ncol(beta)},
#' or a numeric vector of weights with length \code{nrow(beta)}. 
#' @param cell.counts Proportion of cell counts for one cell type in cases
#' where the samples are mainly composed of two cell types (e.g. saliva) (Default: NULL).
#' @param isva Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: TRUE).  
#' @param sva Apply Surrogate Variable Analysis (SVA) to the
#' methylation levels and covariates and include
#' the resulting variables as covariates in a regression model (Default: TRUE).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05). 
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
                        cell.counts=NULL,
                        isva=T, sva=T, ## cate?
                        isva0=F,isva1=F, ## deprecated
                        winsorize.pct=0.05,
                        outlier.iqr.factor=NA, ## typical value = 3
                        most.variable=min(nrow(beta), 50000),
                        featureset=NA,
                        random.seed=20161123,
                        verbose=F) {

    if (isva0 || isva1)
        stop("isva0 and isva1 are deprecated and superceded by isva and sva")
    
    if (is.na(featureset))
        featureset <- guess.featureset(rownames(beta))
    features <- meffil.get.features(featureset)
          
    stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% features$name))
    stopifnot(ncol(beta) == length(variable))
    stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(is.null(weights)
              || is.numeric(weights) && (is.matrix(weights) && nrow(weights) == nrow(beta) && ncol(weights) == ncol(beta)
                                         || is.vector(weights) && length(weights) == nrow(beta)
                                         || is.vector(weights) && length(weights) == ncol(beta)))
    stopifnot(most.variable > 1 && most.variable <= nrow(beta))
    stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && winsorize.pct < 0.5)

    original.variable <- variable
    original.covariates <- covariates
    if (is.character(variable))
        variable <- as.factor(variable)
    
    stopifnot(!is.factor(variable) || is.ordered(variable) || length(levels(variable)) == 2)
    
    msg("Simplifying any categorical variables.", verbose=verbose)
    variable <- simplify.variable(variable)
    if (!is.null(covariates))
        covariates <- do.call(cbind, lapply(covariates, simplify.variable))

    sample.idx <- which(!is.na(variable))
    msg("Removing", ncol(beta) - length(sample.idx), "missing case(s).", verbose=verbose)

    if (is.matrix(weights))
        weights <- weights[,sample.idx]
    if (is.vector(weights) && length(weights) == ncol(beta))
        weights <- weights[sample.idx]
    
    beta <- beta[,sample.idx]
    variable <- variable[sample.idx]

    if (!is.null(covariates))
        covariates <- covariates[sample.idx,,drop=F]

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
    if (!is.null(covariates))
        covariate.sets$all <- covariates

    if (is.numeric(winsorize.pct))  {
        msg(winsorize.pct, "- winsorizing the beta matrix.", verbose=verbose)
        beta <- winsorize(beta, pct=winsorize.pct)
    }

    too.hi <- too.lo <- NULL
    if (is.numeric(outlier.iqr.factor)) {
        q <- rowQuantiles(beta, probs = c(0.25, 0.75), na.rm = T)
        iqr <- q[,2] - q[,1]
        too.hi <- which(beta > q[,2] + outlier.iqr.factor * iqr, arr.ind=T)
        too.lo <- which(beta < q[,1] - outlier.iqr.factor * iqr, arr.ind=T)
        if (nrow(too.hi) > 0) beta[too.hi] <- NA
        if (nrow(too.lo) > 0) beta[too.lo] <- NA
    }

    if (isva || sva) {
        beta.sva <- beta
        
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))
        if (length(autosomal.sites) < most.variable) {
          warning("Probes from the sex chromosomes will be used to calculate surrogate variables.")
        } else {
          beta.sva <- beta.sva[autosomal.sites,]
        }
        var.idx <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable]
        beta.sva <- impute.matrix(beta.sva[var.idx,,drop=F])
        
        if (!is.null(covariates)) 
            mod0 <- model.matrix(~ ., data.frame(covariates, stringsAsFactors=F))
        else
            mod0 <- data.frame(rep(1, length(variable)))
        mod <- cbind(mod0, variable)
        
        if (isva) {
            msg("ISVA.", verbose=verbose)
            set.seed(random.seed)
            isva.ret <- isva(beta.sva, mod, verbose=verbose)
            covariate.sets$isva <- data.frame(covariates, isva.ret$isv, stringsAsFactors=F)
        }
        
        if (sva) {
            msg("SVA.", verbose=verbose)
            set.seed(random.seed)
            sva.ret <- sva(beta.sva, mod=mod, mod0=mod0)
            covariate.sets$sva <- data.frame(covariates, sva.ret$sv, stringsAsFactors=F)
            cat("\n")
        }
    }

    analyses <- sapply(names(covariate.sets), function(name) {
        msg("EWAS for covariate set", name, verbose=verbose)
        covariates <- covariate.sets[[name]]
        ewas(variable,
             beta=beta,
             covariates=covariates,
             batch=batch,
             weights=weights,
             cell.counts=cell.counts,
             winsorize.pct=winsorize.pct)
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
         version=packageVersion("meffil"),
         samples=sample.idx,
         variable=original.variable[sample.idx],
         covariates=original.covariates[sample.idx,,drop=F],
         winsorize.pct=winsorize.pct,
         outlier.iqr.factor=outlier.iqr.factor,
         most.variable=most.variable,
         p.value=p.values,
         coefficient=coefficients,
         analyses=analyses,
         random.seed=random.seed,
         too.hi=too.hi,
         too.lo=too.lo)
}

is.ewas.object <- function(object)
    is.list(object) && "class" %in% names(object) && object$class == "ewas"

#' Test associations between \code{variable} and each row of \code{beta}
#' while adjusting for \code{covariates} (fixed effects) and \code{batch} (random effect).
#' If \code{cell.counts} is not \code{NULL}, then it is assumed that
#' the methylation data is derived from samples with two cell types.
#' \code{cell.counts} should then be a vector of numbers
#' between 0 and 1 of length equal to \code{variable} corresponding
#' to the proportions of cell of a selected cell type in each sample.
#' The regression model is then modified in order to identify
#' associations specifically in the selected cell type (PMID: 24000956).
ewas <- function(variable, beta, covariates=NULL, batch=NULL, weights=NULL, cell.counts=NULL, winsorize.pct=0.05,
                 verbose=F) {
    stopifnot(all(!is.na(variable)))
    stopifnot(length(variable) == ncol(beta))
    stopifnot(is.null(covariates) || nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(is.null(cell.counts)
              || length(cell.counts) == ncol(beta)
              && all(cell.counts >= 0 & cell.counts <= 1))

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

    msg("Linear regression with limma::lmFit", verbose=verbose)
    fit <- NULL
    batch.cor <- NULL
    if (!is.null(batch)) {
        msg("Adjusting for batch effect", verbose=verbose)
        corfit <- duplicateCorrelation(beta, design, block=batch, ndups=1)
        batch.cor <- corfit$consensus
        msg("Linear regression with batch as random effect", verbose=verbose)
        tryCatch(fit <- lmFit(beta, design, block=batch, cor=batch.cor, weights=weights),
                 error=function(e) {
                     print(e)
                     msg("lmFit failed with random effect batch variable, omitting", verbose=verbose)
                 })
    }
    if (is.null(fit)) {
        msg("Linear regression with only fixed effects", verbose=verbose)
        batch <- NULL
        fit <- lmFit(beta, design, weights=weights)
    }
    
    msg("Empirical Bayes", verbose=verbose)
    if (is.numeric(winsorize.pct)) {
      fit.ebayes <- eBayes(fit, robust=T, winsor.tail.p=c(winsorize.pct, winsorize.pct))
    } else { 
      fit.ebayes <- eBayes(fit)
    }

    alpha <- 0.975
    std.error <- (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,"variable"])
    margin.error <- (std.error * qt(alpha, df=fit.ebayes$df.total))
    n <- rowSums(!is.na(beta))

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




