# Shared statistical normalization for BigWig display and future independent R analysis.
# Always supply original unthresholded integer counts, never track values.
capzyme_size_factors <- function(counts, type=c("ratio", "poscounts"),
                                control_rows=NULL, min_features=1L) {
  type <- match.arg(type)
  if (!requireNamespace("DESeq2", quietly=TRUE)) stop("Install DESeq2 for this normalization mode")
  if (!is.matrix(counts) || !is.numeric(counts) || ncol(counts) < 2L || nrow(counts) < 1L ||
      any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)) ||
      any(counts > .Machine$integer.max)) stop("Need nonnegative R-integer-range counts and at least two cohort samples")
  if (any(colSums(counts)==0)) stop("All-zero sample cannot support DESeq2 size-factor estimation")
  if (length(min_features)!=1L || is.na(min_features) || min_features < 1L || min_features!=floor(min_features))
    stop("min_features must be a positive integer")
  rows <- if (is.null(control_rows)) seq_len(nrow(counts)) else control_rows
  if (!is.numeric(rows) || !length(rows) || anyNA(rows) || any(rows!=floor(rows)) ||
      any(rows < 1L | rows > nrow(counts)) || anyDuplicated(rows)) stop("Invalid control row indices")
  chosen <- counts[rows,,drop=FALSE]
  eligible <- if (type=="ratio") rowSums(chosen>0)==ncol(chosen) else rowSums(chosen)>0
  contributions <- colSums(chosen[eligible,,drop=FALSE]>0)
  diagnostics <- list(method=type, requested_features=length(rows),
                      eligible_features=sum(eligible), all_zero_features=sum(rowSums(chosen)==0),
                      zero_entries=sum(chosen==0), contributing_features=contributions)
  if (any(contributions < min_features))
    stop(paste0("Insufficient informative features for ",type," (minimum ",min_features,
                " per sample); contributions: ",paste(names(contributions),contributions,collapse=", "),
                "; all-zero features: ",diagnostics$all_zero_features))
  storage.mode(counts) <- "integer"
  coldata <- data.frame(row.names=colnames(counts), cohort=rep("common",ncol(counts)))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~1)
  if (is.null(control_rows)) {
    dds <- DESeq2::estimateSizeFactors(dds, type=type)
  } else {
    dds <- DESeq2::estimateSizeFactors(dds, type=type, controlGenes=control_rows)
  }
  sf <- DESeq2::sizeFactors(dds)
  if (any(!is.finite(sf) | sf<=0)) stop("DESeq2 returned invalid size factors")
  list(size_factors=sf, diagnostics=diagnostics,
       DESeq2_version=as.character(utils::packageVersion("DESeq2")), R_version=R.version.string)
}
