# Standalone settings, version 2: implemented subset

Use [config/tools.v2.example.yaml](../config/tools.v2.example.yaml) with
`capzyme OPERATION --settings FILE`, optionally `--preview`.
This is a strict, deliberately limited interface. The larger proposed v2 YAML
is not accepted unchanged. Existing direct CLIs and version-1 Snakemake configs
are preserved; `capzyme.py` still accepts version 1 only.

## File and precedence contract

- A sample sheet describes files and metadata; analysis settings describe tools,
  references and requested operation outputs. Keep real sheets/configs private.
- JSON works with Python's standard library. YAML additionally requires PyYAML,
  not Snakemake. Unknown fields, invalid types/choices and duplicate YAML keys
  fail, including in sections not selected by the current operation. Duplicate
  JSON keys are rejected too.
- Explicit CLI options override selected settings, which override documented
  defaults. Boolean options support positive/negative forms, for example
  `--exclude-duplicates` / `--no-exclude-duplicates`. A CLI value of zero is an
  actual override. Unsupported CLI options for an operation are errors.
- Settings-relative paths resolve against the settings file. A sample sheet's
  `input_file` paths resolve against the sheet. Explicit CLI paths resolve against
  the invoking directory. Version-1 working-directory semantics are unchanged.
- Output directories must be new and external. Mapping requires existing indexes;
  neither a missing part nor a hash mismatch triggers an implicit rebuild.
- `--preview` prints resolved absolute paths, effective command, policies and
  ordered metadata without creating results. It can hash/verify a supplied index
  manifest, but is not a substitute for an authorized resource/input preflight.
- Each operation writes `methods.json` with resolved command, precedence, settings
  hash, output, status and adapter hash. Scientific engines retain their own exact
  native commands, versions, reference/input hashes and QC. For BigWigs the record
  also includes the existing engine parser's fully expanded parameters. Trimming
  records Cutadapt's actual executable/version. Preserve these together.

Supported sheet columns: required unique string `sample_id`; `input_file` is
required when selecting reads/BAMs from the sheet. Optional `experiment`,
`treatment`, `biological_source`, `batch`, `display_label`. Row order defines sample
order for counting/default BigWig cohort. Leading zeros in TSV sample IDs remain
intact. Count-only sheets may omit `input_file`. `--sample` selects one library
for trimming/mapping; there is no hidden batch scheduler.

## Support table

| Section/operation | Implemented fields/behavior | Deferred/rejected |
| --- | --- | --- |
| Top level | `schema_version: 2`, optional `experiment_name`, `samples`, `output_dir`, sections below | `execution`, arbitrary extensions |
| `input` | `type: raw_fastq|trimmed_fastq|bam|counts`; `matrix`, `assert_unthresholded` for BigWig | Paired reads, implicit raw→trimmed conversion |
| `references.rna` / `.genome` | `fasta` OR `chrom_sizes` for BigWig; finished FASTA for mapping/counting; `decoys_file` or `no_decoys` for genome mapping | Automatic RNA collapse/viral combination |
| Each `index` | `action: use_existing`, `prefix`, `verification_manifest`; manifest optional; any supplied manifest must match | `build_new`, annotation changes during mapping; use explicit native preparation |
| `trimming` | `adapter`, `quality`, `minimum_length`, `nextseq_trim`; selected only by `trim` | `enabled`, implicit trimming inside mapping |
| `mapping.rna` | `aligner: bowtie2`, `preset: very-sensitive`, `allow_reverse_complement`, `seed`, `minimum_mapq`, `exclude_duplicates`, `threads` | Other aligners/presets, arbitrary extra args |
| `mapping.genome` | `aligner: hisat2|bowtie2`, `preset: very-sensitive`, `seed`, `minimum_mapq`, `exclude_duplicates`, `threads`, `uncertain_cigar: reject|error` | Alternative soft-clipping/orientation policies, arbitrary native args |
| `counting` | `rna_strands`, `genome_strands`: plus/minus/both; `minimum_mapq`, `exclude_duplicates`, `minimum_total_count` | Automatic extra filtered matrix in standalone mode |
| `bigwig.rna` / `.genome` | One `normalization: raw|rpm|deseq2|controls|supplied` per invocation; `strands`, `cohort`, `emit_samples`, `depths`, `matrix_complete_both_strands`, `deseq2_type`, `controls`, `min_control_features`, `size_factors`, `factor_kind`, `factor_metadata` | `normalizations` list, `enabled`, implicit cross-space factors; Implicit BAM selection from a settings sample sheet (use explicit `--bam`/`--bam-manifest`) |
| Index preparation / collapse | Independent native builders or `capzyme collapse`; `index-record` verifies/registers finished index | Settings-driven large reference preparation |
| Snakemake | Existing version-1 routes use the same scientific implementations | Version-2 translation/scheduling (workflow BigWigs use a separate v1 block) |

Valid supported sections for other tools may coexist in a file; only the selected
tool/space consumes its relevant section. An explicitly supplied CLI option for
another operation is rejected. Declarations of input type must match the selected
operation. For example, use `input.type: counts` for browser generation, not the
raw FASTQ document from an earlier stage. `input.matrix` identifies that one
space's matrix; run the other reference space with its own path/settings.

Annotation is a reference-preparation choice: the existing HISAT2 prepare CLI
supports plain GTF **or** splice/exon tables. This subset does not accept those
fields inside `use_existing`; they cannot change embedded annotation of an index.
Native builders and the [index guide](tools/index_preparation.md) expose them.

## Defaults and enforced policies

| Tool | Defaults / enforced behavior |
| --- | --- |
| Trim | Explicit adapter; q20, min15, NextSeq20, one thread |
| RNA alignment | Finished FASTA/existing Bowtie2 index; very-sensitive; reverse complements disabled by default; seed0, MAPQ20, duplicates retained, one thread; flags2308 |
| Genome alignment | HISAT2 default; very-sensitive; no soft clipping; both orientations; seed0/MAPQ20/retain duplicates/one thread; explicit decoys; uncertain CIGAR reject; sort256MiB/worker |
| Count | Both strands unless per-space setting; MAPQ20, flags2308, duplicate/QCfail retention; threshold0 |
| BigWig | Explicit normalization/cohort/unthresholded assertion; both strands; existing BigWig engine defaults (ratio, minimum10 controls, Rscript when needed) |

RNA orientation and count strands are distinct choices. Terminal-coordinate CIGAR
acceptance is a genome-adapter policy, not an invented behavior for external BAMs.
Downstream filtering cannot recover alignments discarded earlier. For an external
BAM choose its own count settings deliberately; metadata is not evidence that
the original alignment used these settings. For complete read accounting use
the verified adapters and preserve original FASTQ/all-record evidence.

```bash
capzyme count --settings "$WORK/bam-tools.yaml" \
  --mapq 30 --no-exclude-duplicates --strand minus \
  --output-dir "$WORK/new-counts" --preview
capzyme bigwig --space genome --settings "$WORK/browser.yaml" \
  --normalization raw --output-dir "$WORK/new-raw" --preview
```

Remove incompatible method-specific settings when changing normalization. For
example `depths` is valid only for RPM; an explicit raw override does not silently
ignore an RPM depth declaration. New files and outputs remain external to Git.

For standalone DESeq2 tracks use [bigwig.deseq2.v2.example.yaml](../config/bigwig.deseq2.v2.example.yaml).
Set `normalization: deseq2`, explicit `deseq2_type: ratio|poscounts`, an ordered
cohort and independent emitted subset. Select an explicit R installation with
CLI `--rscript`; the YAML subset does not accept an Rscript/execution-profile
field. Direct and settings execution share the same method-option checks.
`min_control_features` is rejected outside controls mode; complete-both-strands
assertion applies only to matrix RPM. Existing defaults and valid version-1/direct
scientific behavior are unchanged. See [DESeq2 normalization contract](bigwig2.md).

See [reference/hash contract](reference_verification.md) for observed versus expected identity.

For trimming, explicit CLI numeric/off overrides settings; null settings omit the
NextSeq option, omission defaults to 20, and numeric zero remains enabled.
