# Generate strand-specific BigWigs

**Program:** pyBigWig writes single-base signals; the validated counter prepares
BAM input when requested. Actual R/DESeq2 estimates factors only in selected R
modes. [bigwig2.py](../../workflow/scripts/bigwig2.py) is standalone; the toolbox
exposes matrix and BAM input without changing normalization. Optional version-2
settings primarily describe matrix operations; BAM paths can be supplied explicitly.

```bash
capzyme bigwig --matrix "$INTEGER_MATRIX" \
  --assert-unthresholded --chrom-sizes "$SIZES" --cohort A B \
  --strand both --normalization rpm --depths "$AUDITED_DEPTHS" \
  --output-dir "$WORK/new-rpm-tracks"
capzyme bigwig --space genome \
  --settings "$WORK/browser-tools.yaml" --normalization raw \
  --output-dir "$WORK/new-raw-tracks" --preview

# One command estimates native DESeq2 factors from A+B, rendering only A.
capzyme bigwig --matrix "$INTEGER_MATRIX" \
  --assert-unthresholded --chrom-sizes "$SIZES" --cohort A B --emit-samples A \
  --strand both --normalization deseq2 --deseq2-type ratio \
  --rscript "$RSCRIPT" --output-dir "$WORK/new-deseq2-tracks"
# The same implementation through strict settings (edit paths/cohort first):
capzyme bigwig --settings "$WORK/deseq2.yaml" \
  --rscript "$RSCRIPT" --preview
```

Python/pyBigWig are sufficient for raw/RPM/supplied matrix modes; no R or
Snakemake. BAM input needs SAMtools/BEDTools; deseq2/controls need R/DESeq2.
Inputs are unthresholded nonnegative integer matrices or BAMs, matching unchanged
FASTA/sizes, a declared ordered cohort and selected strands. BAM mode is available through installed `capzyme bigwig --bam FILE --sample-name ID`
(repeat both) or `--bam-manifest FILE`. It uses the same counter and does not
require a matrix assertion. Declare the cohort and matching reference/sizes.

Modes: raw counts; full-precision RPM using audited accepted records; actual
DESeq2 ratio/poscounts; DESeq2 controlGenes; supplied factors/multipliers.
No background subtraction/smoothing/coverage extension/strand negation, automatic
estimator fallback or legacy rounding. Cohort estimation is independent of the
emitted subset. RNA-stage and genome-stage RPM depths are distinct conditional
pools. Never transfer their factors silently. Both strands share a sample's
within-space denominator; it is not a separate plus/minus library size.

Outputs: `tracks/SAMPLE.plus.bw` / `.minus.bw`, `reference.sizes.tsv`,
`normalization_factors.tsv`, `provenance.json`, `commands.log`; R modes add
`normalization/`, BAM mode adds accepted copies/counts. Toolbox adds
`methods.json`. Coordinate p becomes zero-based `[p-1,p)`; strand is relative to
each unchanged FASTA record. Values are nonnegative with documented float32
storage. An otherwise empty track has a declared first-base zero sentinel.

Use the matching custom reference in a genome browser, with a common vertical
range within a space and method. Raw, RPM and DESeq2 units differ. Future R
statistics consume original integers, never tracks. Full normalization, control,
saved-factor and cohort semantics: [BigWig 2 guide](../bigwig2.md).

Start settings from [the DESeq2 example](../../config/bigwig.deseq2.v2.example.yaml).
The estimator must be a deliberate choice; sparse `ratio` failure does not switch
to `poscounts`. RNA and genome require separate invocations/factors. One genomic
factor applies to both strands. DESeq2 values are count / size factor, without
RPM/common-million scaling. Each method preserves a library's spatial shape but
can scale libraries differently. Use shared vertical axes within one method and
space; agreement with native R establishes computation, not biological suitability.
The [normalization contract](../bigwig2.md) describes diagnostics, factor
semantics and limitations.

```bash
capzyme bigwig --bam "$BAM_A" --sample-name A --bam "$BAM_B" --sample-name B \
  --reference "$FASTA" --cohort A B --normalization rpm \
  --output-dir "$WORK/new-bam-rpm"
```

BAM RPM measures its complete accepted both-strand/reference pool internally.
Workflow generation instead reuses the existing unthresholded matrices; see
[optional workflow tracks](../workflow_bigwigs.md).
