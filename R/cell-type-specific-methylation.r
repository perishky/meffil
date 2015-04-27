meffil.cell.type.specific.methylation <- function(beta, cell.types, num.sites=50, verbose=F) {
    msg("cell types", paste(unique(cell.types), collapse=", "), verbose=verbose)

    stopifnot(is.matrix(beta))
    stopifnot(ncol(beta) == length(cell.types))
    stopifnot(num.sites > 0 && num.sites < nrow(beta))

    number.cell.types <- length(unique(cell.types))
    design <- model.matrix(~ 0 + cell.types)

    msg("fitting linear model with limma::lmFit", verbose=verbose)
    fit <- limma::lmFit(beta, design=design)
    contrasts <- diag(x=2, nrow=number.cell.types) - 1
    contrasts[which(contrasts == 1)] <- nrow(contrasts)-1

    msg("computing t-statistics with limma::eBayes", verbose=verbose)
    fit.eb <- limma::eBayes(contrasts.fit(fit, contrasts))
    
    sites.idx <- apply(fit.eb$t, 2, function(t.statistics) {        
        c(order(t.statistics, decreasing=F)[1:num.sites],
          order(t.statistics, decreasing=T)[1:num.sites])
    })
    sites.idx <- unique(as.vector(sites.idx))

    ## checked that this wwas working as expected:
    ## 1. should be equal 
    ## 6*fit$coefficients[1000,1] - sum(fit$coefficients[1000,-1])
    ## fit.eb$coefficients[1000,1]

    ## 2. should be large (methylation differences between cell type and the rest)
    ## quantile(apply(abs(fit.eb$coefficients[sites.idx,])/6, 1, max, na.rm=T))

    ## 3. p-values should all be really small
    ## quantile(apply(fit.eb$p.value[sites.idx,], 1, min, na.rm=T))

    ## 4. t-statistics should be positive or negative and all far from 0.
    ## quantile(apply(fit.eb$t[sites.idx,], 1, min, na.rm=T))
    ## quantile(apply(fit.eb$t[sites.idx,], 1, max, na.rm=T))
    ## quantile(apply(abs(fit.eb$t[sites.idx,]), 1, max, na.rm=T))

    ret <- fit$coefficients[sites.idx,]
    colnames(ret) <- gsub("cell.types", "", colnames(ret))
    ret
}

