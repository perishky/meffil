#' Concordance between genotypes and SNP betas
#'
#' genotypes <- meffil.extract.genotypes(dataset.prefix)
#' snp.betas <- meffil.snp.betas(qc.objects)
#' meffil.snp.concordance(snp.betas, genotypes[rownames(snp.betas),colnames(snp.betas)])
#'
#' @return Returns a list of two vectors:
#' one providing concordances between genotypes and SNP betas for matched samples,
#' a second providing concordances between genotypes and SNP betas for matched SNPs.
#' 
#' @export
meffil.snp.concordance <- function(snp.betas, genotypes, snp.threshold=0.99) {
    stopifnot(length(colnames(snp.betas)) == ncol(snp.betas))
    stopifnot(length(rownames(snp.betas)) == nrow(snp.betas))
    stopifnot(all(colnames(snp.betas) == colnames(genotypes)))
    stopifnot(all(rownames(snp.betas) == rownames(genotypes)))
 
    beta.genotypes <- calculate.beta.genotypes(snp.betas)

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
    snp.idx <- which(snp.concordance > snp.threshold)

    if (length(snp.idx) > 0)
        sample.concordance <- colSums(beta.genotypes[snp.idx,,drop=F] == genotypes[snp.idx,,drop=F],
                                      na.rm=T)/length(snp.idx)
    else
        sample.concordance <- rep(NA, ncol(genotypes))

    list(sample=sample.concordance,
         snp=snp.concordance)
}

calculate.beta.genotypes <- function(snp.betas, centers=c(0.2,0.5,0.8)) {
    t(apply(snp.betas,1,function(x) kmeans(x, centers=centers)$cluster - 1))
}


