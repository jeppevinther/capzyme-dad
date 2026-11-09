# Choose a tool

Start with the [standalone tool catalogue](toolbox.md): explicit inputs, outputs,
upstream programs and copyable commands for each operation. The optional
[version-2 settings adapter](tool_settings.md) works without Snakemake.

The remainder of this guide covers **optional version-1 Snakemake coordination**.
[Getting started](getting_started.md) provides its synthetic RNA-only tutorial,
installation notes and working-directory conventions. Standalone tools remain
the primary interfaces; coordination uses their shared scientific implementations.

In examples, `CAPZYME_REPO` is the absolute checkout path, `WORK` is an external
project directory and all uppercase input variables refer to files you supply.
Create configurations in `$WORK`, use absolute data paths, and activate the
appropriate software environment before running. Public templates contain `null`
placeholders: fill them in a text editor. Keep private configuration, research
data and results outside Git.

## Choose a route

| Starting material / objective | Route | Config template | Orchestrated by Snakemake? |
| --- | --- | --- | --- |
| Raw single-end FASTQ → RNA-reference counts | `rna` (recommended start) | [core](../config/core.example.yaml) or [fixed](../config/fixed.example.yaml) | Yes |
| Raw FASTQ → RNA reference, then genome for aligner-unmapped reads | `sequential` | RNA config plus `genome` block from [genome template](../config/genome.example.yaml) | Yes |
| Already trimmed FASTQ → genome counts | `genome` | [genome](../config/genome.example.yaml) | Yes |
| Existing single-end RNA-reference BAMs → counts | `bam` | [existing BAM](../config/existing_bam.example.yaml) | Yes; standalone counter also available |
| Counts/BAMs → BigWigs | `capzyme bigwig` or optional workflow tracks | CLI, optional v2 settings, or v1 `bigwigs` block | Optional |

`check`, `preview` and `run` use the same `--route` and configuration. All launcher
routes invoke existing validators/rules. Direct Snakemake remains supported:
`workflow/Snakefile` serves RNA/sequential; `workflow/genome.smk` serves genome-only;
`workflow/existing_bam.smk` serves BAM counting. The launcher prints its exact
Snakemake command. Use at least four cores for the original RNA streaming mapping
pipeline; see [resources](resources.md) for budgeting and enforcement limits.

## Raw FASTQ to RNA-reference counts

Required: raw single-end Phred33 FASTQ per library; TSV `sample_id`, `fastq`;
RNA reference FASTA; explicit adapter/trimming parameters; new external outputs.
The RNA mapping default is Bowtie2 `--very-sensitive` with `norc: true` and seed
zero. `norc: true` disallows reverse-complement alignments relative to each FASTA
record. It does not infer biological library strandedness.

Generated mode collapses ncRNAs sharing their configured 5′ prefix (default
30 nt), writes membership/representative FASTA, combines optional ordered viral
records and builds a managed Bowtie2 index. Review the
[reference strategy](reference_strategy.md). Fixed mode uses an existing FASTA
and complete six-part `.bt2` or `.bt2l` index. It automatically records hashes, enforces optional supplied pins, and
verifies reference/index content before reads are processed. It does not regenerate
the source reference or index; mixed suffix families fail.

```bash
capzyme workflow check --route rna --config "$WORK/rna.yaml"
capzyme workflow preview --route rna --config "$WORK/rna.yaml"
capzyme workflow run --route rna --config "$WORK/rna.yaml" \
  --cores 4 --memory-mb 16384
```

Main results: `trimmed/`, Cutadapt `qc/`, accepted `bam/`, managed `reference/`,
`counts/5ends.raw.tsv`, optional `counts/5ends.min_total_N.tsv`, logs and provenance.
“Raw” means unthresholded integer counts **after alignment filtering**. It does
not mean unfiltered alignments or raw-read coverage.

## Sequential RNA-reference → genome mapping

Required: the raw FASTQ/RNA configuration above, plus an existing complete genome
index, matching FASTA, explicit decoy policy and a `genome` config block. Choose
HISAT2 for splice-containing genomes or Bowtie2 for unspliced genomes. The recommended
Conda environment includes both aligners and their builders/inspectors. Relevant pins are documented
in [sequential.yaml](../workflow/envs/sequential.yaml); fresh deployment of that
per-rule environment is not claimed as tested here.

```bash
capzyme workflow doctor --route sequential --aligner hisat2
capzyme workflow check --route sequential --config "$WORK/sequential.yaml"
capzyme workflow preview --route sequential --config "$WORK/sequential.yaml" \
  --cores 8 --memory-mb 65536
capzyme workflow run --route sequential --config "$WORK/sequential.yaml" \
  --cores 8 --memory-mb 65536
```

Only reads the **RNA aligner classified unmapped** are forwarded. A low-MAPQ RNA
alignment is still a mapping; forwarding it as though unmapped would allow an
ambiguous read to be assigned twice or change the reference-selection policy.
Mapped rejects remain in the RNA all-record archive. Handoff IDs, sequences and
qualities are verified at each stage; final accounting now verifies their complete
provenance/checksum connection and the BAM/QC artifacts it consumes.

HISAT2 retains `--very-sensitive --no-softclip` and both orientations. Genomic
acceptance also excludes declared decoys, any S/H/P, and terminal I/D/N; internal
I/D and splice N remain allowed. This ensures the terminal aligned nucleotide
can define a conservative 5′ coordinate. Mapping and generic counting policies
are detailed in [genome mapping](genome_mapping.md).

Main results: each `rna/SAMPLE/` and `genome/SAMPLE/` has accepted/all-record BAMs,
unmapped FASTQ, QC, dictionary, command log and provenance/publication receipt.
`qc/SAMPLE.sequential.json` records exclusive raw/trimmed/stage totals, zero accepted
read-ID overlap and verified handoff hashes. Keep `counts/5ends.raw.tsv` (RNA)
and `counts/genome.5ends.raw.tsv` (genome) separate. Their references represent
different spaces; identical coordinate numbers do not identify the same feature.
No normalization factors are transferred across them. See
[sequential details](sequential_mapping.md) for resources/restarts.

## Genome-only, already trimmed reads

Required: already trimmed single-end FASTQs; TSV `sample_id`, `fastq`; matching
FASTA and existing complete genome index; explicit decoy policy. The config's
`genome` block is the same as above. This route performs no trimming, RNA mapping
or index construction.

```bash
capzyme workflow check --route genome --config "$WORK/genome.yaml"
capzyme workflow preview --route genome --config "$WORK/genome.yaml" \
  --cores 8 --memory-mb 65536
capzyme workflow run --route genome --config "$WORK/genome.yaml" \
  --cores 8 --memory-mb 65536
```

Outputs are `genome/SAMPLE/`, `counts/genome.5ends.raw.tsv` and
`provenance/genome.json`. Standalone `bowtie2_genome.py map` and
`hisat2_genome.py map` remain available and accept no Snakemake configuration.
Their `prepare` operations create separate managed indexes. HISAT2 preparation
can use compatible GTF or pre-extracted splice/exon tables; verify assembly,
identifier, strand and coordinate compatibility first. Preparing a human annotated
index needs a separate resource plan and is not part of mapping or the tutorial.
See [preparation commands](genome_mapping.md#hisat2-genome).

Decoys must be designated explicitly: supply `genome.decoys` with one exact
identifier per line, or declare `no_decoys: true` for a genuinely decoy-free
reference. HISAT2 requires that choice; the workflow requires it for either
genome aligner. Keep decoys in the original reference/index dictionary, excluding
their alignments through the declared acceptance policy rather than guessing
names or altering the index.

## Count existing BAMs

Required: one or more single-end BAMs and the matching reference FASTA. The BAM
workflow manifest is `sample_id`, `bam`, with desired column order. The template
contains an optional historical comparison: **omit `comparison` or use `{}` for
ordinary counting**; do not leave the historical target's `null` fields in place.
Supply your own FASTA hash. Sorting/filtering status is validated; previously
discarded reads cannot be restored.

```bash
capzyme workflow doctor --route bam
capzyme workflow check --route bam --config "$WORK/bam.yaml"
capzyme workflow preview --route bam --config "$WORK/bam.yaml"
capzyme workflow run --route bam --config "$WORK/bam.yaml"

# Independent FivePrime Count 2, two ordered BAMs, no Snakemake:
capzyme count --reference "$FASTA" \
  --bam "$BAM_A" --sample-name A --bam "$BAM_B" --sample-name B \
  --strand both --output-dir "$WORK/new-counts"
```

Both use the same BEDTools counting engine. Default: MAPQ ≥20, exclude flags 2308,
retain duplicates/QC-failure records, no row-total threshold. Standalone counting
uses new output/log paths. Sources are not resorted/reindexed; filtering uses
copies and sorting is skipped for actually ordered copies. Workflow results are
accepted `bam/`, `reference/sizes.tsv`, raw/optional thresholded `counts/`, logs and
provenance. See [counter details](existing_bam_counting.md).

## Standalone BigWig visualization

BigWig generation is optional in all workflow routes through the `bigwigs`
block; see [workflow BigWigs](workflow_bigwigs.md). Alternatively use the standalone tool
with unthresholded integer counts or accepted BAMs and the matching reference/
sizes dictionary. It additionally requires pyBigWig 0.3.26. R/DESeq2 is needed
only for `deseq2` and `controls`; environment specifications and full CLI examples
are in [BigWig 2](bigwig2.md). These are optional routes after the beginner tutorial.

```bash
capzyme bigwig --help
capzyme bigwig \
  --matrix "$RAW_COUNTS" --assert-unthresholded --chrom-sizes "$SIZES" \
  --cohort A B --normalization raw --output-dir "$WORK/new-raw-tracks"
```

| Mode | Meaning / extra inputs |
| --- | --- |
| `raw` | Original integer 5′ counts, stored as track values |
| `rpm` | Counts × 1,000,000 / accepted alignment depth; use verified depths if the matrix is not the complete both-strand universe |
| `deseq2` | Actual R/DESeq2 size factors from the declared full sample cohort; choose `ratio` or `poscounts` explicitly when appropriate |
| `controls` | Actual DESeq2 `controlGenes` estimation using a TSV of stable 5′ positions; factors apply to all positions |
| `supplied` | Previously calculated sample-specific size factors/multipliers; saved metadata supports verified track regeneration |

Outputs include separate plus/minus `.bw` tracks, `normalization_factors.tsv`,
provenance and logs. BigWig coordinates are zero-based `[p-1,p)` for matrix
position p. Strands are relative to the supplied reference sequence, including
reverse-orientation viral records. Tracks are single-nucleotide ends without
smoothing or read-coverage extension. DESeq2-normalized values are not RPM.
Future independent statistical R/DESeq2 analysis consumes the original integer
matrix, never BigWig values, and does not depend on BigWig generation.

All routes currently support single-end data, not paired-fragment counting,
multiline FASTQ, automatic UMI deduplication or a complete differential-expression
analysis. For completion,
logs, recovery and failures, use [results/troubleshooting](results_and_troubleshooting.md).

Optional browser tracks are described in [workflow BigWigs](workflow_bigwigs.md).

See [reference/hash contract](reference_verification.md) for observed versus expected identity.
