#' Functional normalization
#'
#' Apply functional normalization to a set of Infinium HumanMethylation450 BeadChip IDAT files.
#'
#' Fortin JP, Labbe A, Lemire M, Zanke BW, Hudson TJ, Fertig EJ, Greenwood CM, Hansen KD.
#' Functional normalization of 450k methylation array data
#' improves replication in large cancer studies.
#' Genome Biol. 2014 Dec 3;15(12):503. doi: 10.1186/s13059-014-0503-2.
#' PMID: 25599564
#'
#' @param samplesheet Output from \code{\link{meffil.read.samplesheet}()} or
#' \code{\link{meffil.create.samplesheet}()}.
#'
#' Arguments to \code{\link{meffil.qc}()}:
#' @param number.quantiles
#' @param detection.threshold
#' @param bead.threshold
#' @param sex.cutoff
#' @param featureset
#' @param chip
#' @param cell.type.reference
#'
#' Argument to \code{\link{meffil.qc.summary}()}:
#' @param qc.parameters (parameters)
#'
#' Arguments to \code{\link{meffil.qc.report}()}:
#' @param qc.file (output.file)
#' @param author
#' @param study
#'
#' Arguments to \code{\link{meffil.normalize.quantiles}()}:
#' @param number.pcs
#' @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
#' along with control matrix principal components (Default: NULL).
#' @param random.effects Names of columns in samplesheet that should be included as random effects
#' (Default: NULL).
#'
#' Arguments to \code{\link{meffil.normalize.samples}()}:
#' @param just.beta
#' @param pseudo
#'
#' Arguments to \code{\link{meffil.normalization.summary}()}:
#' @param norm.parameters (parameters)
#' @param norm.file (output.file)
#'
#' Other:
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#'
#' @return A list:
#' - qc.summary \code{\link{meffil.qc.summary}()} output.
#' - norm \code{\link{meffil.normalize.quantiles}()} output.
#' - beta Normalized beta matrix (methylation levels) of class \code{\link[bigmemory]{big.matrix}}.
#' - norm.summary \code{\link{meffil.normalization.summary}()} output.
#'
#' @export
meffil.normalize.dataset <- function(samplesheet,
                                     
                                     ## meffil.qc
                                     number.quantiles=500,
                                     detection.threshold=0.01,
                                     bead.threshold=3,
                                     sex.cutoff=-2,
                                     featureset=NULL,
                                     chip=NULL,
                                     cell.type.reference=NULL,
                                     
                                     ## meffil.qc.summary
                                     qc.parameters=meffil.qc.parameters(),
                                     
                                     ## meffil.qc.report
                                     qc.file = "meffil-qc-report.md",
                                     author = "Analyst",
                                     study = "IlluminaHuman450 data",

                                     ## meffil.normalize.quantiles
                                     number.pcs=2,
                                     fixed.effects=NULL,
                                     random.effects=NULL,

                                     ## meffil.normalize.samples
                                     pseudo=100,
                                     just.beta=T,

                                     ## meffil.normalization.summary
                                     norm.parameters=NULL,
                                     norm.file="meffil-normalization-report.md",
                                     
                                     verbose=FALSE) {
    
    qc.objects <- meffil.qc(samplesheet,
                            number.quantiles=number.quantiles,
                            verbose=verbose,
                            detection.threshold=detection.threshold,
                            bead.threshold=bead.threshold,
                            sex.cutoff=sex.cutoff,
                            featureset=featureset,
                            chip=chip,
                            cell.type.reference=cell.type.reference)
    
    qc.summary <- meffil.qc.summary(qc.objects,
                                    parameters=qc.parameters,
                                    verbose=verbose)


    meffil.qc.report(qc.summary,
                     output.file=qc.file,
                     author=author,
                     study=study)

    if (nrow(qc.summary$bad.samples) > 0)
        qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
    
    norm.objects <- meffil.normalize.quantiles(qc.objects,
                                               fixed.effects=fixed.effects,
                                               random.effects=random.effects,
                                               number.pcs=number.pcs,
                                               verbose=verbose)
    
    norm <- meffil.normalize.samples(norm.objects,
                                     pseudo=pseudo,
                                     just.beta=just.beta,
                                     cpglist.remove=qc.summary$bad.cpgs$name,
                                     verbose=verbose)
    if (just.beta)
        beta <- norm
    else
        beta <- get.beta(norm$M, norm$U)
    
    if (is.null(norm.parameters))
        norm.parameters <- meffil.normalization.parameters(norm.objects)

    norm.summary <- meffil.normalization.summary(beta,
                                                 norm.objects=norm.objects,
                                                 parameters=norm.parameters,
                                                 verbose=verbose)

    meffil.normalization.report(norm.summary,
                                output.file=norm.file,
                                author=author,
                                study=study)

    list(qc.summary=qc.summary,
         norm.objects=norm.objects,
         M=if (just.beta) NULL else norm$M,
         U=if (just.beta) NULL else norm$U,
         beta=beta,
         norm.summary=norm.summary)
}    


