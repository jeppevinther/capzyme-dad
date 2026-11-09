# Sequential RNA-reference and genome mapping

Add a `genome` mapping to the established FASTQ configuration to enable the
optional second stage. Without it, the original RNA-only rules and matrix paths
remain active. The [genome template](../config/genome.example.yaml) documents
the complete optional block; copy that block into a generated/fixed RNA-reference
configuration. References, indexes, policy files, raw data and outputs are external.
The manifest remains `sample_id`, `fastq` in the desired matrix-column order.

```bash
nice -n 10 snakemake --snakefile workflow/Snakefile \
  --configfile "$PRIVATE_SEQUENTIAL_CONFIG" --executor local \
  --cores 32 --resources mem_mb=98304 --printshellcmds
```

CPU and memory values are execution settings, not hardcoded workflow-wide limits.
Use a detached session and external process-tree/host monitoring. `mem_mb` is a
Snakemake scheduling budget, not an operating-system memory cap. Launch-time
checks must cover available RAM, load, storage and applicable quota. Preserve a
host available-memory reserve of max(64 GiB, 10% of total RAM). Mapping rules
reserve one thread for the I/O helper; each standalone tool aligns and sorts
sequentially with the remaining CPU allocation. Counting/accounting request two
threads for their parent and subprocess. Numerical library threads should be one.

## Data flow and scientific policy

1. Cutadapt produces the established trimmed FASTQ and JSON report.
2. `rna_handoff.py` runs Bowtie2 with the established very-sensitive preset,
   configured `--norc`, seed, MAPQ and SAM-flag settings. Its accepted BAM uses
   the same SAMtools filtering as the original route. It retains an all-record
   BAM and the Bowtie2 `--un` FASTQ before any MAPQ/flag filtering.
3. A disk-backed read-ID/fingerprint check requires that handoff records are
   exactly the aligner's unmapped inputs, preserving sequence and quality.
   Mapped reads rejected by MAPQ or flags stay outside the handoff.
4. The genome rule invokes the existing `bowtie2_genome.py map` or
   `hisat2_genome.py map` command. HISAT2 retains its fixed
   `--very-sensitive --no-softclip` policy. Genome acceptance retains its MAPQ,
   flags-2308, explicit decoy and conservative CIGAR safeguards. Duplicates and
   QC-failure reads remain eligible unless duplicate exclusion is requested.
5. `sequential_accounting` joins primary identities and checks disjoint accepted
   read sets, exhaustive handoff and exclusive outcomes across both stages. It
   requires complete/input-verified stage provenance, validates publication
   receipts and consumed BAM/QC hashes, and requires recorded RNA handoff-output,
   recorded genome-input and current handoff checksums to agree. Consumed hashes
   are rechecked before a successful report can be written. Matching IDs alone
   do not establish sequence/quality integrity.
6. Both spaces use the unchanged `core.matrix` / BEDTools counter, the same
   implementation used by standalone FivePrime Count 2. RNA and genomic counts
   are separate original integer matrices with reference-relative strands and
   one-based `start == end`. Both unthresholded matrices are retained.

For the human route, choose the separately validated annotated index and explicit
audited decoy list. All six/eight index parts, matching FASTA, decoy file and
`genome.policy_files` are declared DAG inputs. The mapper fingerprints its actual
inputs before/after execution; cohort provenance fingerprints policy files too.
Index generation is deliberately external to this mapping DAG.

## Outputs and restart behavior

| Artifact | Location under external `output_dir` |
| --- | --- |
| RNA accepted/all-record BAMs | `rna/SAMPLE/{accepted,alignments}.bam` |
| Exact first-stage handoff | `rna/SAMPLE/unmapped.fastq` |
| Genome accepted/all-record BAMs | `genome/SAMPLE/{accepted,alignments}.bam` |
| Final aligner-unmapped reads | `genome/SAMPLE/unmapped.fastq` |
| Per-stage QC, commands and provenance | Each stage directory: `qc.json`, `commands.jsonl`, `provenance.json`, workflow `publication.json` |
| Dictionaries and genome per-reference QC | `reference.sizes.tsv`, genome `mapping_by_reference.tsv` |
| Exclusive identity/accounting check | `qc/SAMPLE.sequential.json` |
| RNA integer matrix | `counts/5ends.raw.tsv` (optional existing thresholded output too) |
| Genomic integer matrix | `counts/genome.5ends.raw.tsv` |
| Cohort provenance | `provenance/{run,genome}.json` |

The optional route always retains all-record BAMs; `keep_diagnostic_sam` belongs
to the original RNA-only route. It does not control these BAMs. No genomic BAM
index is required for counting, so none is generated automatically.

Each tool invocation uses a fresh output in `.attempts/SPACE/SAMPLE/`. Only a
successful, input-verified result is published. Failed attempts remain available
for diagnosis. A publication receipt links the original standalone output path
to its final workflow location and pins the unmodified tool provenance. Previous
workflow-owned results are moved to `.history/` when a
stage is rerun; unmanaged nonempty directories are refused. Snakemake may remove
declared output files before a forced rerun, so `.history/` is a recovery aid,
not a full versioned backup of every prior file. Preserve external run directories
when immutable complete run archives are needed. No automatic cleanup is applied.

Use `--rerun-incomplete` after diagnosing/correcting a failed run. A normal
completed restart should report `Nothing to be done`; changing an index part or
declared policy input invalidates dependent mapping and counting rules.

## Independent entry points

Genome-only Snakemake orchestration accepts already trimmed FASTQs and the
[genome-only configuration](../config/genome.example.yaml):

```bash
nice -n 10 snakemake --snakefile workflow/genome.smk \
  --configfile "$PRIVATE_GENOME_CONFIG" --executor local \
  --cores 32 --resources mem_mb=98304
```

It skips trimming, RNA mapping and reference construction. The two existing
[standalone genome commands](genome_mapping.md) remain independently usable.
The RNA handoff adapter is also independently executable:

```bash
python workflow/scripts/rna_handoff.py --help
python workflow/scripts/rna_handoff.py --fastq "$TRIMMED_FASTQ" \
  --reference "$RNA_FASTA" --index "$RNA_INDEX" --norc \
  --mapq 20 --threads 8 --output-dir "$NEW_RNA_OUTPUT"
```

Normalization and differential-expression analysis are separate downstream
operations. Factors are not transferred between these reference dictionaries.

For a guided starting route, use [Getting started](getting_started.md). The
[practical guide](usage.md) compares all routes, and
[results/troubleshooting](results_and_troubleshooting.md) explains checksum errors,
normal resume and history limits. Final accounting verifies the full handoff
sequence/quality and artifact checksum chain; matching read IDs alone is insufficient.
