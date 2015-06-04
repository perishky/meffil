#' Extract genotype data from PLINK .raw files for Illumina 450K SNPs
#' 
#' @param filenames A vector of filenames of PLINK .raw files from which to extract
#' genotype data.
#' @return Matrix with rows corresponding to SNPs, columns to samples and values
#' equal to 0, 1 or 2 corresponding to genotypes.
#' @examples
#' R> writeLines(meffil.snp.probes(), con="snp-names.txt")
#' shell> plink --bfile dataset --extract snp-names.txt --recodeA --out genotypes.raw --noweb
#' R> filenames <- "genotypes.raw"
#' R> genotypes <- meffil.extract.genotypes(filenames)
#' 
#' @export
meffil.extract.genotypes <- function(filenames, verbose=F) {
    stopifnot(all(sapply(filenames, file.exists)))
    
    table.list <- lapply(filenames, function(filename) {
        msg("Reading plink file", filename, verbose=verbose)
        read.table(filename, header=T)
    })

    sample.names <- table.list[[1]][,1] ## family id
    genotypes <- lapply(table.list, function(genotype.table) {
        genotype.table[match(sample.names, genotype.table[,1]), ## match family ids
                       -(1:6),
                       drop=F]
    })
    genotypes <- do.call(cbind, genotypes)
    colnames(genotypes) <- gsub("_[GCAT]{1}$", "", colnames(genotypes)) 
    rownames(genotypes) <- sample.names

    t(as.matrix(genotypes))
}




