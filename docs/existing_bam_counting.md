# Standalone and existing-BAM 5′-end counting

The standalone counter accepts one or multiple existing single-end BAMs without
remapping. See the [quick command guide](tools/counting.md) and the exact
comparison procedure below. Keep original BAMs unchanged.

The standalone counter and both Snakemake counting routes use `core.matrix` and
`core.threshold`. BEDTools `genomecov -5 -dz` defines the aligned 5′ position;
there is no second CIGAR or scientific counting implementation. Important
operations should remain independently executable, with Snakemake coordinating
the same implementation.

## Standalone command

Requires Python 3, SAMtools and BEDTools on PATH. No Python packages beyond the
standard library are needed for the standalone commands. Run `--help` for options.

```bash
python "$REPO/workflow/scripts/count_5ends.py" --help
python "$REPO/workflow/scripts/count_5ends.py" \
  --reference "$REFERENCE" --reference-sha256 "$REFERENCE_SHA256" \
  --sample-name B --sample-name A \
  --strand both --mapq 20 --min-total-count 10 \
  --output "$NEW_OUTPUT/counts.tsv" "$BAM_B" "$BAM_A"
```

One or more BAMs are positional arguments. Repeat `--sample-name` in the desired
order, once per BAM; without names, BAM filename stems are used. Names must be
unique and may contain letters, numbers, dots, underscores and hyphens (starting
with a letter or number). Duplicate source paths and reserved coordinate names
are rejected. This is a single-end counter, not a paired-fragment counter.

Defaults: MAPQ >=20, exclude flags 2308 (unmapped, secondary, supplementary),
retain duplicate and QC-fail flags, both strands, no total-count threshold.
`--exclude-duplicates` also excludes flag 1024. `--strand plus|minus|both` and
`--min-total-count N` control strand selection and the inclusive row total across
all supplied samples. No multimapping weight or additional deduplication is used.

The output TSV has `chromosome`, `strand`, `start`, `end`, then the ordered sample
columns. Coordinates are one-based; `start == end`. Counts are raw integers and
absent sample/feature combinations are zero. Rows have deterministic natural
reference order, position order and strand order. Text order may differ from the
historical script while keyed scientific values remain identical.

All inputs are read-only. The reference dictionary must match each BAM. The tool
checks BAM structure, reapplies filtering into an output-local copy and streams
the filtered alignments to validate actual coordinate order against the BAM
dictionary. Already ordered copies skip sorting, regardless of their header;
only unordered copies are sorted. This avoids redundant sorting/indexing and
trusting stale sort headers or unspecified prior filters. The validation adds
a sequential decoding pass; it cannot restore previously removed alignments.
When needed, sorting uses one thread and 256 MiB buffers. No source-adjacent FASTA
or BAM indexes are written. Temporary files are removed when the operation exits.

The output and its JSON-lines log must be new paths distinct from every input.
The default log is `OUTPUT.jsonl`; `--log` selects another new path. It records
resolved tool versions, command arguments, command exit status/wall time, input
sizes/timestamps, reference/output checksums and counting settings. BAM content
hashes are not computed automatically. Source size and modification time are
checked before and after preparation; synthetic tests additionally check full
source byte hashes. This is not a cryptographic guarantee against concurrent
source replacement; arrange stable inputs for real runs.

## Snakemake entry point

Use `workflow/existing_bam.smk`, with the portable template
[existing_bam.example.yaml](../config/existing_bam.example.yaml). The original
`workflow/Snakefile`, its FASTQ configuration and rules are unchanged.
`REPO` denotes this checkout; the other example variables below denote external
input/output locations selected in private configuration.

The BAM manifest is a tab-separated file with exactly `sample_id` and `bam`.
Manifest row order becomes matrix column order. Paths resolve relative to the
Snakemake working directory; absolute paths are recommended in private server
configurations. Reference SHA-256 expectations are optional; observations are recorded automatically. No Bowtie2 index is required.
Input/reference/research data must be outside the checkout and the output tree.

The DAG makes output-local reference sizes, prepares each BAM using the same
`prepare_bam` function as the standalone tool, assembles raw counts using
`core.matrix`, applies `core.threshold`, and writes provenance. The optional
`comparison` section runs an exact historical comparison. It contains no trimming,
reference collapse, mapping or index-building jobs. BAM preparations can run in
parallel; matrix assembly runs once. Concurrency and total memory scheduling
limits come from the command line; per-rule requests come from the config.

```bash
# Synthetic data are authorized; real inputs require separate execution approval.
snakemake --snakefile "$REPO/workflow/existing_bam.smk" \
  --configfile "$PRIVATE_CONFIG" --directory "$EXTERNAL_WORKDIR" \
  --executor local --cores 4 --resources mem_mb=16384 --dry-run
```

For an authorized real-data run after checking resources and arranging active
memory/load monitoring, remove `--dry-run` and use the local execution policy:

```bash
nice -n 10 snakemake --snakefile "$REPO/workflow/existing_bam.smk" \
  --configfile "$PRIVATE_CONFIG" --directory "$EXTERNAL_WORKDIR" \
  --executor local --cores 4 --resources mem_mb=16384 --printshellcmds
```

Under the external `output_dir`, the workflow retains:

- `bam/SAMPLE.accepted.bam`: filtered, coordinate-ordered copies; originals unchanged.
- `reference/sizes.tsv`: derived reference sizes; no source-side index.
- `counts/5ends.raw.tsv` and `counts/5ends.min_total_N.tsv`: raw union and
  inclusive total-threshold matrices, with sample columns in manifest order.
- `logs/` and `provenance/run.json`: commands, statuses, parameters and provenance.
- `comparison/exact.json`: optional exact comparison report. A failed comparison
  also retains its report in `logs/compare.jsonl`.

The public example contains the RSV comparison section. For ordinary counting
without a historical target, omit `comparison` entirely (or use `{}`); supply
your own reference checksum and input paths. `resources.cores`/`resources.mem_mb`
record intended budgets in configuration; the command-line `--cores` and
`--resources mem_mb=...` establish the actual Snakemake scheduling limits.

With `--use-conda`, the coordinator still needs Snakemake 9.9.0 and its supported
Python environment. The new preparation environment pins Python 3.12/SAMtools
1.22; the counting environment pins BEDTools 2.31.1. The completed historical
server comparison used the verified installed SAMtools 1.21/BEDTools 2.31.0 on
PATH without `--use-conda`; they were reverified before execution. Do not silently
substitute tools between comparison arms.

## Exact matrix comparison

The comparator is independently executable and is also used by Snakemake:

```bash
python "$REPO/workflow/scripts/compare_matrices.py" \
  --expected "$HISTORICAL_MATRIX" --observed "$NEW_MATRIX" \
  --sizes "$OUTPUT_REFERENCE_SIZES" --min-total-count 10 \
  --expected-features 11580 --expected-sha256 "$HISTORICAL_MATRIX_SHA256" \
  --report "$NEW_OUTPUT/comparison.json"
```

Keys are reference identifier, strand and one-based position; `start == end` is
validated. Sample IDs and their order must match exactly. Counts must be
nonnegative integers, keys unique, coordinates inside the dictionary and row
totals above the inclusive threshold. Compare feature sets, every integer count
and sample totals at zero tolerance. Byte equality is reported separately from
keyed equality. A mismatch exits nonzero and records aggregate differences plus
bounded private examples. Malformed matrices also fail rather than being fixed.
The Snakemake comparison log retains the report even when Snakemake removes a
failed output. Neither command reorders sample columns, repairs coordinates or
changes threshold to obtain agreement.

See [reference/hash contract](reference_verification.md) for observed versus expected identity.
