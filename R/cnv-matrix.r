
#' Create matrix of CNV values
#'
#' @param cnv Output from \code{\link{meffil.calculate.cnv}()}.
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}.
#' @return Matrix of ncpg x nsample
#' @export
meffil.cnv.matrix <- function(cnv, featureset="450k") {
    features <- meffil:::cnv.features(featureset)
    features$id <- with(features, paste(as.character(chromosome), position))

    cnv.matrix <- sapply(cnv, function(segments)  {
        idx <- order(factor(segments$chrom, levels=unique(features$chrom)), segments$loc.start)
        segments <- segments[idx,]
        segments$start.id <- with(segments, paste(as.character(chrom), as.integer(loc.start)))
        segments$end.id <- with(segments, paste(as.character(chrom), as.integer(loc.end)))
        
        start.idx <- match(segments$start.id, features$id)
        end.idx <- match(segments$end.id, features$id)
        good.idx <- which(!is.na(start.idx) & !is.na(end.idx))
        start.idx <- start.idx[good.idx]
        end.idx <- end.idx[good.idx]
        feature.idx <- unlist(lapply(1:length(start.idx), function(i) start.idx[i]:end.idx[i]))
        
        scores <- rep(NA, nrow(features))
        names(scores) <- features$name
        scores[feature.idx] <- with(segments[good.idx,], rep(as.numeric(seg.mean),
                                                             end.idx-start.idx+1))
        scores
    })
    cnv.matrix
}
           
