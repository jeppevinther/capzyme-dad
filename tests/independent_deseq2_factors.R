# Verification only: direct native DESeq2 on the original coordinate matrix.
# Does not source the CapZyme normalization helper or implement an estimator.
a <- commandArgs(trailingOnly=TRUE)
if (length(a)<4L) stop('Usage: Rscript --vanilla independent_deseq2_factors.R MATRIX NEW_OUTPUT_DIR SAMPLE1 SAMPLE2 ...')
if (!requireNamespace('DESeq2', quietly=TRUE)) stop('Select an Rscript environment containing DESeq2')
if (dir.exists(a[2]) || file.exists(a[2])) stop('Choose a new output directory')
dir.create(a[2], recursive=TRUE)
options(digits=17)
suppressPackageStartupMessages(library(DESeq2))
capture.output(sessionInfo(),file=file.path(a[2],'session.txt'))
cohort <- a[-c(1L,2L)]
t <- read.delim(a[1], check.names=FALSE, quote='', comment.char='')
if (!identical(colnames(t),c('chromosome','strand','start','end',cohort))) stop('Original matrix/header/cohort order mismatch')
x <- as.matrix(t[,cohort,drop=FALSE]); rm(t); invisible(gc())
if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x<0) ||
    any(x!=floor(x)) || any(x>.Machine$integer.max)) stop('Invalid original integer counts')
# Diagnostics are independent eligibility observations, not an estimator.
eligible <- rowSums(x>0)==ncol(x)
contributing <- colSums(x[eligible,,drop=FALSE]>0)
diagnostics <- data.frame(sample_id=cohort, method='ratio', requested_features=nrow(x),
                          eligible_features=sum(eligible), contributing_features=unname(contributing),
                          all_zero_features=sum(rowSums(x)==0), zero_entries=sum(x==0),
                          accepted_records=unname(colSums(x)),
                          R_version=as.character(getRversion()), DESeq2_version=as.character(packageVersion('DESeq2')))
write.table(diagnostics,file.path(a[2],'diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if (any(contributing<1L) || any(colSums(x)==0)) stop('ratio has insufficient informative features; no fallback estimator was selected')
d <- DESeqDataSetFromMatrix(x,data.frame(row.names=cohort,cohort=rep('all',length(cohort))),~1)
d <- estimateSizeFactors(d,type='ratio')
sf <- sizeFactors(d)
if (!identical(names(sf),cohort) || any(!is.finite(sf) | sf<=0)) stop('Invalid direct DESeq2 size factors')
write.table(data.frame(sample_id=names(sf),size_factor=unname(sf)),file.path(a[2],'factors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
