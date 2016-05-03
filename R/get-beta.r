#' Infinium HumanMethylation450 BeadChip methylation levels
#'
#' Compute beta values (methylation levels) from methylated/unmethylated signals
#'
#' @param M Methylated signal matrix.
#' @param U Unmethylated signal matrix.
#' @param pseudo Value to add to the denominator to make the methylation estimate more stable.
#' @return Matrix of 0..1 methylation level estimates.
#' Equal to methylated/(methylated + unmethylated + pseudo).
#'
#' @export
meffil.get.beta <- function(M, U, pseudo=100) {
    stopifnot(is.matrix(M)
              && is.matrix(U)
              && nrow(M) == nrow(U)
              && ncol(M) == ncol(U))
    get.beta(M,U,pseudo)
}

get.beta <- function(M,U,pseudo=100) {
    M[]/(M[]+U[]+pseudo)
}
