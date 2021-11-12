add.surrogate.variables <- function(beta, sites, samples,
                                    most.variable,
                                    variable, covariates,
                                    winsorize.pct,
                                    outlier.iqr.factor,
                                    random.seed, isva, sva, smartsva, n.sv, smartsva.alpha=0.25, verbose=F) {

    stopifnot(length(samples) == length(variable))
    stopifnot(is.null(covariates) || nrow(covariates) == length(samples)) 
    
    meffil:::msg("Calculating CpG variance", verbose = verbose)
    var.sites <- meffil.most.variable.cpgs(beta, n=most.variable, sites=sites, samples=samples, winsorize.pct=winsorize.pct, outlier.iqr.factor=outlier.iqr.factor)
  
    if (is.matrix(beta)) {
        beta <- beta[var.sites,,drop=F]
        if (!is.null(samples))
            beta <- beta[,samples,drop=F]
    }
    else {
        beta <- retrieve.gds.methylation(beta, var.sites, samples) 
    }
    beta <- meffil.handle.outliers(beta, winsorize.pct, outlier.iqr.factor)
    beta <- meffil:::impute.matrix(beta, margin=1)    

    if (!is.null(covariates)) {
        cov.frame <- model.frame(~., data.frame(covariates, stringsAsFactors=F), na.action=na.pass)
        mod0 <- model.matrix(~., cov.frame)
    }
    else
        mod0 <- matrix(1, ncol=1, nrow=length(variable))
    mod <- cbind(mod0, variable)

    covariate.sets <- list()
    last.ret <- NULL
    if (isva) {
        meffil:::msg("ISVA.", verbose=verbose)
        set.seed(random.seed)
        isva.ret <- isva(beta, mod, ncomp=n.sv, verbose=verbose)
        if (!is.null(covariates))
            covariate.sets$isva <- data.frame(covariates, isva.ret$isv, stringsAsFactors=F)
        else
            covariate.sets$isva <- as.data.frame(isva.ret$isv)
        last.ret <- isva.ret
        cat("\n")
    }
        
    if (sva) {
        meffil:::msg("SVA.", verbose=verbose)
        set.seed(random.seed)
        sva.ret <- sva(beta, mod=mod, mod0=mod0, n.sv=n.sv)
        if (!is.null(covariates))
            covariate.sets$sva <- data.frame(covariates, sva.ret$sv, stringsAsFactors=F)
        else
            covariate.sets$sva <- as.data.frame(sva.ret$sv)
        last.ret <- sva.ret
        cat("\n")
    }

    if (smartsva) {
        meffil:::msg("smartSVA.", verbose=verbose)
        set.seed(random.seed)
        if (is.null(n.sv)) {
            beta.res <- t(resid(lm(t(beta) ~ ., data=as.data.frame(mod))))
            n.sv <- EstDimRMT(beta.res, FALSE)$dim + 1
        }
        smartsva.ret <- smartsva.cpp(beta, mod=mod, mod0=mod0, n.sv=n.sv, alpha=smartsva.alpha)                
        if (!is.null(covariates))
            covariate.sets$smartsva <- data.frame(covariates, smartsva.ret$sv, stringsAsFactors=F)
        else
            covariate.sets$smartsva <- as.data.frame(smartsva.ret$sv)
        last.ret <- smartsva.ret
        cat("\n")
    }
    
    list(sets=covariate.sets, last=last.ret)
}
