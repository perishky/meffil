#' Infinium HumanMethylation450 BeadChip methylation levels
#'
#' Compute beta values (methylation levels) from methylated/unmethylated signals
#'
#' @param mu Methylated/unmethylated signal from \code{\link{meffil.rg.to.mu}()}
#' or \code{\link{meffil.normalize.sample}()}.
#' @param pseudo Value to add to the denominator to make the methylation estimate more stable.
#' @return Vector of 0..1 methylation level estimates.
#' Equal to methylated/(methylated + unmethylated + pseudo).
#'
#' @export
meffil.get.beta <- function(mu, pseudo=100) {
    mu$M/(mu$M+mu$U+pseudo)
}
