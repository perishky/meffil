#' Concordance between genotypes and SNP betas
#'
#' Calculating cordance between methylation-derived genotypes and user-provided genotypes 
#' is a bit complicated because some beadchip probes might not perform well and some samples 
#' provide poor-quality data. We require, when calculating concordances, 
#' that only high-quality probes for high-quality samples are included in calculations. 
#' Formally, a high-quality SNP is a SNP for which `snp.threshold`% of the genotypes (in high-quality samples) 
#' determined by DNA methylation match those provided by the user.  
#' A high-quality sample is a sample for which `sample.threshold`% of the genotypes (from high-quality probes) 
#' determined by DNA methylation match those provided by the user.  
#'
#' @param snp.betas Matrix of methylation-derived genotypes (rows=SNPs, columns=samples)
#' @param genotypes Matrix of user-provided genotypes (rows=SNPs, columns=samples)
#' @param snp.threshold Minimum percentage (of high-quality) samples whose user-provided genotypes for a SNP must 
#' match the methylation-derived genotypes to be called high-quality (Default: 0.99). 
#' @param sample.threshold Minimum percentage (of high-quality) SNPs whose user-provided genotypes for a sample must
#' match the methylation-derived genotypes to be called high-quality (Default: 0.9).
#' @return Returns a list of two vectors:
#'     - one providing concordances between genotypes and SNP betas for matched samples,
#'     - a second providing concordances between genotypes and SNP betas for matched SNPs.
#' as well as the genotype matrix derived from 'snp.betas'.
#' @examples
#' genotypes <- meffil.extract.genotypes(raw.filenames)
#' snp.betas <- meffil.snp.betas(qc.objects)
#' meffil.snp.concordance(snp.betas, genotypes[rownames(snp.betas),colnames(snp.betas)])
#'
#' @export
meffil.snp.concordance <- function(snp.betas, genotypes,
                                   snp.threshold=0.99,
                                   sample.threshold=0.9) {
    stopifnot(length(colnames(snp.betas)) == ncol(snp.betas))
    stopifnot(length(rownames(snp.betas)) == nrow(snp.betas))
    stopifnot(all(colnames(snp.betas) == colnames(genotypes)))
    stopifnot(all(rownames(snp.betas) == rownames(genotypes)))
 
    beta.genotypes <- calculate.beta.genotypes(snp.betas)

    ## calculate snp concordances
    snp.concordance <- sapply(rownames(genotypes), function(snp) {
        beta.factor <- factor(beta.genotypes[snp,],levels=0:2)
        geno.factor <- factor(genotypes[snp,],levels=0:2)
        counts <- table(beta.factor, geno.factor)
        diag1 <- sum(counts[1,1] + counts[2,2] + counts[3,3])
        diag2 <- sum(counts[3,1] + counts[2,2] + counts[1,3])
        if (diag1 < diag2) {
            beta.genotypes[snp,] <<- 2 - beta.genotypes[snp,]
            diag2/sum(counts)
        }
        else
            diag1/sum(counts)
    })

    col.concordance <- function(x,y,rows=rep(T,nrow(x)))
        colMeans((x==y)[which(rows),,drop=F], na.rm=T)
    row.concordance <- function(x,y,cols=rep(T,ncol(x)))
        rowMeans((x==y)[,which(cols),drop=F], na.rm=T)
    n.true <- function(x) sum(x, na.rm=T)
    
    ## calc sample concordances
    sample.concordance <- col.concordance(beta.genotypes, genotypes)

    if (ncol(genotypes) >= 20) {
        ## assume that the top 50% of samples are not mismatched
        good.smp <- sample.concordance >= median(sample.concordance,na.rm=T)
        old.smp <- good.smp
                
        for (i in 1:20) {
            snp.concordance <- row.concordance(beta.genotypes, genotypes, good.smp)
            good.snp <- snp.concordance >= snp.threshold
            sample.concordance <- col.concordance(beta.genotypes, genotypes, good.snp)
            good.smp <- sample.concordance >= sample.threshold

            if (n.true(good.snp) < 2 || n.true(good.smp) < 2) {
                snp.concordance <- row.concordance(beta.genotypes, genotypes)
                sample.concordance <- col.concordance(beta.genotypes, genotypes)
                break
            }

            if (all(good.smp == old.smp, na.rm=T))
                break

            old.smp <- good.smp
        }
    }
    
    list(sample=sample.concordance,
         snp=snp.concordance,
         beta.genotypes=beta.genotypes)
}

calculate.beta.genotypes <- function(snp.betas, centers=c(0.2,0.5,0.8)) {
    x <- t(apply(snp.betas,1,function(x) {
        tryCatch(kmeans(x, centers=centers)$cluster - 1,
                 error=function(e) {
                     cluster <- rep(1,ncol(snp.betas))
                     cluster[which(x < min(centers))] <- 0
                     cluster[which(x > max(centers))] <- 2
                     cluster
                 })
    }))
    dimnames(x) <- dimnames(snp.betas)
    x
}


