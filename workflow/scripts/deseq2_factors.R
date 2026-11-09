# Internal CLI adapter; the reusable R function itself has no BigWig dependency.
a <- commandArgs(trailingOnly=TRUE)
if (length(a)!=7L) stop("Expected helper, count TSV, control-index file or -, method, minimum, output TSV, diagnostics TSV")
source(a[1])
x <- read.delim(a[2], row.names=1, check.names=FALSE, quote="", comment.char="")
x <- as.matrix(x)
controls <- if (a[3]=="-") NULL else scan(a[3], what=integer(), quiet=TRUE)
answer <- capzyme_size_factors(x, type=a[4], control_rows=controls, min_features=as.integer(a[5]))
options(digits=17)
write.table(data.frame(sample_id=names(answer$size_factors),size_factor=unname(answer$size_factors)),
            a[6], sep="\t", quote=FALSE,row.names=FALSE)
d <- answer$diagnostics
write.table(data.frame(sample_id=names(answer$size_factors), method=d$method,
            contributing_features=unname(d$contributing_features),
            eligible_features=d$eligible_features, requested_features=d$requested_features,
            all_zero_features=d$all_zero_features, zero_entries=d$zero_entries,
            DESeq2_version=answer$DESeq2_version, R_version=answer$R_version),
            a[7], sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(), file=paste0(a[7],".session.txt"))
