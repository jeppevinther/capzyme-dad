# Results and troubleshooting

Start with the [beginner tutorial](getting_started.md); use [practical usage](usage.md)
to choose a route. Keep configuration, input data, results and Snakemake metadata
in external project storage. In commands below, `CONFIG` is your external config
and `CAPZYME_REPO` is the absolute source checkout. Add the same `--route` used
for your run (`rna` is the default).

## Recognise completion and find the evidence

Successful execution exits zero and Snakemake reports all scheduled steps done.
Inspect the matrix and QC/provenance too; the existence of a BAM alone does not
prove the entire workflow completed. A completed normal restart reports
`Nothing to be done`. A dry run only describes work; it has not validated real
alignment/count output.

| Route | Principal files under `output_dir` |
| --- | --- |
| RNA-only | `bam/SAMPLE.accepted.bam` and `.bai`; `trimmed/`; `reference/`; `qc/SAMPLE.cutadapt.json`; `counts/5ends.raw.tsv`; `provenance/run.json` |
| Sequential | `rna/SAMPLE/` and `genome/SAMPLE/`; `qc/SAMPLE.sequential.json`; RNA `counts/5ends.raw.tsv`; `counts/genome.5ends.raw.tsv`; `provenance/run.json` and `genome.json` |
| Genome-only | `genome/SAMPLE/`; `counts/genome.5ends.raw.tsv`; `provenance/genome.json` |
| Existing BAM | Output-local `bam/SAMPLE.accepted.bam`; `reference/sizes.tsv`; raw/thresholded `counts/`; `provenance/run.json` |
| Standalone BigWig | `tracks/SAMPLE.plus.bw` / `.minus.bw`; `normalization_factors.tsv`; `provenance.json`; `commands.log` |

Workflow logs live in `logs/`. Stage directories in sequential/genome mapping
retain `alignments.bam` (all records), `accepted.bam`, `unmapped.fastq`, `qc.json`,
`commands.jsonl`, `provenance.json` and workflow `publication.json`. Tool command
logs can contain JSON interleaved with native tool messages; do not assume every
line is JSON. Snakemake's coordinator log is in the external workdir's
`.snakemake/log/`. Existing-BAM/standalone counter logs record filtering, sorting
decisions and counting. Keep full logs external: they can contain private paths.

The count TSV starts with `chromosome`, `strand`, `start`, `end`, followed by
sample columns in manifest order. Values are nonnegative integers; `start=end`
is one-based. Missing sample/feature combinations have zero counts. A thresholded
matrix intentionally omits low-total features; the raw/unthresholded matrix is
retained and is the input for normalization/statistical analysis. BigWigs use
zero-based intervals and floating visualization values.

## Reconcile reads

Cutadapt JSON records input and retained reads. For RNA-only runs, combine that
with Bowtie2 logs and accepted BAM counts; a single final report covering all
read outcomes/MultiQC is not yet implemented for that route. The synthetic tutorial
has 18 raw reads, 17 trimmed/accepted, and matrix totals B=12, A=5, empty=0.

Sequential `qc/SAMPLE.sequential.json` provides:

- `raw_input_reads`, `trimmed_input_reads` and `trimming_removed_reads`.
- `exclusive_raw_outcomes`: trimming removals plus RNA accepted/rejected and
  genome accepted/rejected/final-unmapped outcomes, summing to raw input reads.
- `forwarded_unmapped_reads`: a **subtotal**, equal to genome-stage input depth;
  do not add it to exclusive outcomes a second time.
- `accepted_read_id_overlap: 0`, required for the two accepted sets.
- `handoff_integrity`: recorded RNA output, recorded genome input and current
  handoff SHA-256; all three must agree before verification is reported true.
- `verified_artifacts_sha256`: the consumed stage provenance/publication, BAM and
  QC artifacts, checked before and after identity accounting.

The stages' primary-read outcomes form an exclusive partition. Secondary and
supplementary alignment **records** are reported separately and do not represent
additional input reads. An unmapped read has no aligner mapping and can pass to
the next reference. A mapped low-MAPQ, decoy, flag or uncertain-CIGAR rejection
is still a mapped read; it is not silently recycled as unmapped. Duplicates and
QC-failure reads remain eligible under the established default flags.

Final verification now rejects missing/incomplete provenance, checksum mismatch,
or inconsistent publication. The receipt links a published stage to its original
fresh attempt without modifying standalone provenance. Matching IDs alone cannot
establish sequence/quality identity. Hash receipts protect against mixed or
changed artifacts; they are not digital signatures authenticating an untrusted
third party. Earlier pilot reports predate this corrected final checksum gate;
see [the sequential integrity contract](sequential_mapping.md).

## Resume after interruption

First preserve the workdir and outputs and read the last error. Restore the
intended environment, stable inputs and configuration, then repeat the same
`run` command and budgets. Finished, current jobs are reused.

```bash
capzyme workflow run --config "$CONFIG" --cores 4 --memory-mb 16384
# If Snakemake reports incomplete files after an interruption:
capzyme workflow run --config "$CONFIG" --cores 4 --memory-mb 16384 \
  --rerun-incomplete
```

Keep `--route sequential|genome|bam` and any `--workdir` override when applicable.
Do not regenerate the synthetic inputs or delete `.snakemake` as a routine way to
resume. If a lock is reported, first verify that no workflow is still running;
unlocking is an expert recovery action, not an automatic retry. The launcher
does not unlock or force jobs.

## Diagnose common failures

| Symptom | What to check / recover |
| --- | --- |
| `python`, `snakemake` or another tool missing | Activate the intended isolated environment; run `doctor` for your route. Read its executable paths. Genome/BigWig/R dependencies are additional to the RNA tutorial environment. |
| Invalid YAML, disallowed `null`, unknown field or manifest error | Use the right route/template; fill required fields. TSV headers and sample IDs must be exact and unique; paths with spaces should be quoted in shell commands. `check` performs existing schema validation. |
| Missing file / unexpected relative path | Confirm the launcher's printed workdir. Default is the config directory; direct Snakemake differs unless `--directory` is supplied. Prefer absolute input paths. Inputs cannot live under `output_dir` or inside the checkout. |
| Reference checksum/dictionary or index-family error | Select the intended assembly/FASTA and all six/eight matching index parts. Do not mix small/large families, rename incompatible contigs, regenerate over a historical index or edit checksum pins to bypass an unexplained mismatch. Prepare a separate index if required. |
| Handoff or artifact integrity error | Stop. Preserve logs and original files. Restore a mutually matching, trustworthy complete stage set, or rerun from verified inputs into new outputs. Changing only a FASTQ, BAM, QC or provenance hash cannot repair the chain. |
| Standalone output already exists | Choose a new output directory (genome/BigWig) or new matrix/log paths (counter). Standalone tools deliberately refuse replacement. Snakemake's managed attempts handle workflow restarts. |
| Killed process, out of memory, disk/quota error, heavy host load | Stop and inspect host/process resources. Free-space totals may not reveal project quota. Increase approved budgets/reservations or reduce concurrency after diagnosis; `nice` is not a memory cap. Do not automatically switch shared servers. |
| DESeq2 ratio fails on sparse data | Read diagnostics; there may be no informative positive-across-all-sample features. Choose `poscounts` explicitly only if appropriate; no automatic estimator substitution occurs. Controls require sufficiently informative stable positions. |

Treat edits to input files, tool versions, policy files, code or configuration as
potentially changing the analysis. Preview the new DAG. Use a **new config/output
directory** for a changed scientific analysis or for a complete immutable record.
Preserve the original raw inputs and old results. Resource-only changes can also
change scheduling/recorded parameters; inspect the preview rather than assuming
every change triggers or avoids reruns. No automatic cleanup of attempts/history
or old session files is performed.

## Fresh attempts, publication and history

Mapping tools require a fresh output directory. Workflow adapters run them under
`.attempts/SPACE/SAMPLE/attempt-.../result`, then publish only a complete, verified
result. Standalone provenance retains the original attempt output path, and
`publication.json` records the destination and provenance hash. Failed attempts
remain for diagnosis. A successful published result moves out of its attempt.

On rerun, remaining previous workflow-owned stage content moves into `.history/`.
Unmanaged nonempty destinations are refused. **History is not a complete backup**:
Snakemake may remove declared outputs before invoking the adapter, so an older
BAM/FASTQ may already be gone when residual files enter history. Preserve a whole
external run separately when you need a complete immutable archive. Avoid manually
moving individual stage files or editing their publication records.

## Advanced force reruns

Use normal resume first. Direct Snakemake `--forcerun RULE` deliberately reruns a
selected rule; `--forceall` can rerun the entire selected DAG. These options can
remove/recreate declared output files, regenerate managed references in generated
mode, publish new attempts and invalidate downstream counts. Retained `.history`
does not guarantee recovery of removed outputs. Fixed mode has no reference-build
rules, but forcing its analysis jobs can still replace generated analysis outputs.

The launcher intentionally offers no force flag. Experts should copy the exact
printed Snakemake command, review a dry run with the intended force option, and
obtain any required experimental execution approval before running it. Do not
use forcing, altered filters or edited metadata to manufacture agreement with an
old result.

Optional browser tracks are described in [workflow BigWigs](workflow_bigwigs.md).
