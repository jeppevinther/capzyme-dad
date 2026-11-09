`deseq2_normalization.R` exposes the reusable `capzyme_size_factors(counts, type,
control_rows, min_features)` function. It uses actual DESeq2 on original integer
counts (coordinates removed, sample columns preserved). `type` is ratio or
poscounts; optional one-based control row indices use controlGenes. The result
contains size factors, sparse-data diagnostics and R/DESeq2 versions.

BigWig 2's internal adapter sources this same function. Independent R analysis
can source it without Python/BigWig or a display-normalization prerequisite.
No differential-expression pipeline or models are implemented yet. Never use
normalized tracks as the R count input. See docs/bigwig2.md for the API contract.

For standalone tracks use explicit `--normalization deseq2 --deseq2-type ratio`
and `--rscript` with the complete ordered factor cohort. RNA and genomic spaces
use separate estimates; genomic strands share a sample factor. Track units are
integer counts / native DESeq2 size factor, without additional rescaling.

`tests/independent_deseq2_factors.R` is a verification fixture that calls native
DESeq2 directly on an original coordinate matrix without sourcing this helper.
It records ratio eligibility diagnostics before estimation. It is not a second
production normalization component. The [normalization guide](../docs/bigwig2.md)
distinguishes computational diagnostics from biological suitability. Differential
testing remains a separate future component.
