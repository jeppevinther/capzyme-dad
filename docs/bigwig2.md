<!-- Direct script syntax remains supported for compatibility. -->
# BigWig 2 usage

For a short standalone introduction, see [Choose a tool](toolbox.md) and the
[BigWig quick guide](tools/bigwigs.md). Optional version-2 **matrix-input** settings
are implemented in `capzyme bigwig`; the direct CLI below remains unchanged,
including BAM input and all five normalization modes. No Snakemake is needed.

`workflow/scripts/bigwig2.py` is a standalone single-end 5′-end visualization tool.
It shares unchanged BAM preparation and BEDTools counting code with FivePrime
Count 2. There is no background subtraction, smoothing, coverage extension,
strand negation, `legacy-rpm` mode or common-million scale.

Use Python 3.12 and pyBigWig 0.3.26. BAM inputs additionally need SAMtools and
BEDTools. Only `deseq2`/`controls` invoke Rscript and actual DESeq2; simple and
supplied modes require no R. The alternative environments
`workflow/envs/bigwig.yaml` and `bigwig_deseq2.yaml` document these dependencies.
`--rscript` can select an external R installation. No package installation is
performed by the command itself. See `--help` for all options.

## Inputs and normalization

Supply exactly one input route: `--matrix`, repeated `--bam` with ordered
`--sample-name`, or `--bam-manifest` with exactly `sample_id` and `bam`. Supply
one matching `--reference` FASTA or two-column `--chrom-sizes` TSV. Reference
identifiers must be unique, lengths positive and positions valid. BAM dictionaries
must match. FASTA/sizes and original BAMs remain unchanged; preparation uses
copies, sorts only unordered copies and builds no source indexes.

Matrix columns are `chromosome`, `strand`, one-based `start`, `end`, then samples.
Counts must be original nonnegative integers, keys unique and `start=end`.
`--assert-unthresholded` is required: the caller certifies that no count features
were thresholded away. This cannot be inferred from numeric contents alone. Use
FivePrime Count 2's unthresholded output or Snakemake `counts/5ends.raw.tsv`, not
the historical threshold-10 matrix. Never pass normalized BigWig-derived values.

Declare an ordered `--cohort` for estimation; `--emit-samples` only chooses
outputs, without changing factors. All references and both alignment strands in
the estimation matrix share each sample factor. BAM input is counted on both
strands internally even if only one strand is emitted. Matrix input uses the
whole supplied feature universe, which must be chosen consistently across samples.

| Method | Values | Important details |
| --- | --- | --- |
| `raw` | Original counts | Empty libraries can emit empty tracks |
| `rpm` | count × 1,000,000 / accepted alignment records | Full precision before BigWig float storage; one depth across strands/contigs |
| `deseq2` | count / DESeq2 size factor | Actual R estimator; units are normalized counts, **not RPM** |
| `controls` | count / control-based DESeq2 size factor | Same factors apply to every position, not only controls |
| `supplied` | count / supplied size factor, or count × supplied multiplier | Explicit `--factor-kind`; no refitting/re-centering |

RPM counts primary mapped records passing MAPQ >=20 and `-F 2308`; duplicate and
QC-fail records remain. `--mapq` and `--exclude-duplicates` can change the declared
BAM policy. Each record contributes an aligned 5′ end, not a paired fragment or a
unique molecule. MAPQ 255 passes. Matrix RPM needs
`--matrix-complete-both-strands`, certifying the full-reference accepted-end union,
or `--depths`: `sample_id`, positive integer `accepted_records`, common `scope_id`.
A plus-only/cropped matrix sum is not implicitly the historical library depth.
Zero depth is an error. Sequential ncRNA→genome stage depths are conditional;
use an explicitly audited disjoint common pool if external depths are supplied.
The tool cannot establish disjointness from an arbitrary depth TSV.

DESeq2 modes support `--deseq2-type ratio|poscounts` (default ratio), at least two
cohort samples, and nonempty libraries. The R helper calls
`DESeq2::estimateSizeFactors`, with `controlGenes` in controls mode. No pseudocounts,
Python estimator, automatic method switch, dispersion fitting or differential
expression test are used. The R function accepts counts within R integer range.
It reports eligible/informative features, zero entries, contributions per sample,
estimator type and R/DESeq2 versions. Sparse ratio failures are explicit; choose
poscounts deliberately if biologically appropriate.
[DESeq2 documentation](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

Controls TSV has exactly `chromosome`, `strand`, `start`, `end`; no annotation is
needed. Invalid bounds/strands, nonsingle-base keys and duplicate keys fail.
Valid controls absent from the observed union are reported as missing and
excluded as uninformative; observed all-zero controls are reported separately.
R requires `--min-control-features` informative controls **per sample** (default
10, configurable but at least 2). This computational safeguard cannot prove
biological stability of selected monophosphate positions. All-zero samples,
insufficient controls or invalid factors fail with diagnostics, without fallback.

## Examples

`REPO` is the checkout; other variables denote external inputs/new output paths.
These examples do not authorize experimental execution.

```bash
python "$REPO/workflow/scripts/bigwig2.py" --help

# Raw matrix: neither BAM processing nor R needed.
python "$REPO/workflow/scripts/bigwig2.py" \
  --matrix "$RAW_MATRIX" --assert-unthresholded --chrom-sizes "$SIZES" \
  --cohort Rpp NudC MDE untreated --normalization raw --output-dir "$NEW_RAW"

# BAM RPM; one library is valid. Its denominator includes both strands.
python "$REPO/workflow/scripts/bigwig2.py" \
  --bam "$BAM" --sample-name Rpp --reference "$FASTA" \
  --cohort Rpp --normalization rpm --strands plus --output-dir "$NEW_RPM"

# One command estimates actual DESeq2 factors and produces tracks.
python "$REPO/workflow/scripts/bigwig2.py" \
  --matrix "$RAW_MATRIX" --assert-unthresholded --chrom-sizes "$SIZES" \
  --cohort Rpp NudC MDE untreated --normalization deseq2 \
  --deseq2-type ratio --rscript "$RSCRIPT" --output-dir "$NEW_DESEQ"

# Control-based estimation; applies factors to all positions.
python "$REPO/workflow/scripts/bigwig2.py" \
  --matrix "$RAW_MATRIX" --assert-unthresholded --chrom-sizes "$SIZES" \
  --cohort Rpp NudC MDE untreated --normalization controls \
  --controls "$MONOPHOSPHATE_CONTROLS" --deseq2-type ratio \
  --min-control-features 10 --output-dir "$NEW_CONTROLS"

# Reuse an entire saved factors table while regenerating only MDE from one BAM.
# The full saved cohort is retained; missing other sample BAMs are not needed.
python "$REPO/workflow/scripts/bigwig2.py" \
  --bam "$MDE_BAM" --sample-name MDE --reference "$FASTA" \
  --normalization supplied --size-factors "$NEW_DESEQ/normalization_factors.tsv" \
  --factor-kind size-factor --factor-metadata "$NEW_DESEQ/provenance.json" \
  --emit-samples MDE --output-dir "$NEW_MDE"
```

## Outputs and reuse

The output directory must be new and outside the source checkout. It contains:

- `tracks/SAMPLE.plus.bw` / `.minus.bw`, selected by `--strands` (default both).
- `normalization_factors.tsv`: sample ID, size factor where relevant, effective
  multiplier, depth where relevant, method, units and cohort digest. `supplied`
  accepts this file or an external TSV with `sample_id` and `size_factor` or
  `multiplier`. All factors must be finite and positive; rows unique. Reuse the
  matching `--factor-kind` to avoid accidental double normalization.
- `provenance.json`: input/output/code hashes, reference dictionary, declared
  cohort, emitted samples, parameters/units, source count digests, preparation,
  normalization diagnostics and versions. Failed runs retain failure provenance.
- `commands.log`: tool commands, statuses and full R failure messages.
- `normalization/` in R modes: raw cohort count table, control indices when used,
  original R factors, diagnostics and R session info. BAM mode also retains
  accepted copies and `counts.raw.tsv`; experimental outputs belong on the server.

A saved provenance JSON verifies the factor-file checksum, cohort, reference
names/lengths and emitted samples' raw count digests, permitting single-sample
regeneration with unchanged data. Bare external factors are accepted with an
explicit declared cohort (or factor-table order) and a warning that source
context is not verified. Metadata verification does not authenticate an arbitrary
third-party export; preserve trustworthy original artifacts.

Position p becomes `[p-1,p)`. Both tracks are nonnegative; a plus alignment on a
reverse-orientation viral reference stays plus on **that same reference**. No
cross-contig coordinate conversion occurs. BigWig stores floating values, so
float32 rounding is expected. Sparse absent positions are missing coverage. An
empty strand has a zero sentinel at the first contig's first base, explicitly
listed in provenance. Nonzero values outside representable float32 range fail.
Use the matching custom RNA/viral assembly in IGV or UCSC.
[pyBigWig writer API](https://github.com/deeptools/pyBigWig).

## Reusable R component and orchestration boundary

`R/deseq2_normalization.R` exposes `capzyme_size_factors(raw_counts, type,
control_rows, min_features)`. Independent R analysis can source it directly,
using the original integer matrix with coordinate columns removed and correct
sample order. It has no Python, BigWig or display dependency. BigWig's internal
R adapter calls the same function. No differential-expression workflow is
implemented here, and statistical analyses never consume BigWig values.

The optional [workflow rules](workflow_bigwigs.md) now call the standalone CLI; cohort factor estimation must precede per-sample
expansion, or a single invocation can emit all tracks. The existing FASTQ and
BAM-to-count workflows and their scientific functions remain unchanged.

## Standalone DESeq2 interface checks

The direct CLI and optional `capzyme bigwig` use the same option validator and
normalization implementation. Both accept explicit `--deseq2-type`, ordered
`--cohort`, independent `--emit-samples`, strands/dictionary and `--rscript`.
[Settings example](../config/bigwig.deseq2.v2.example.yaml): select the R executable
with CLI `--rscript`; an `rscript` YAML field is not part of the current subset.
No Snakemake or separately run normalization script is required.

Incompatible estimator/depth/factor/control options fail before processing.
`--min-control-features` is controls-only; its default remains 10.
`--matrix-complete-both-strands` is matrix-RPM-only. An Rscript environment locator
may still be specified in simple modes without invoking R, preserving the direct
interface's existing R-independent use. Missing executable errors identify
`--rscript`; failed R commands retain their dependency/feature errors in
`commands.log`, with no estimator fallback. Ordinary DESeq2 estimation retains
the existing minimum of one informative feature per sample, with a warning below
ten. This computational safeguard is not a biological suitability criterion.

Actual R/DESeq2 diagnostics and session information are retained in
`normalization/r_diagnostics.tsv` and its `.session.txt` companion. Source matrix
hashes, sample count digests, cohort/emit order, estimator, exact commands, factor
table/track hashes and units remain in provenance. DESeq2 and RPM preserve each
library's within-track spatial shape but can produce different relative scaling
between libraries; their absolute units differ. No re-centering or common-million
scaling is performed after native DESeq2 factor estimation.

See the [BigWig command guide](tools/bigwigs.md) for installed inputs and examples.
Retain native-R diagnostics and independently decoded interval checks for each validation.
