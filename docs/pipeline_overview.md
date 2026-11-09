# Pipeline architecture

CapZyme-seq exposes independent tools through `capzyme`; optional Snakemake
workflows coordinate those same implementations. Standalone commands do not
require Snakemake. Start with the [tool guide](toolbox.md) and [installation](installation.md).

## Independent operations

| Command | Input | Output / shared implementation |
| --- | --- | --- |
| `capzyme collapse` | RNA FASTA | Collapsed FASTA and membership table; preserved prefix-collapse implementation |
| Native `bowtie2-build` / `hisat2-build` | FASTA, optional compatible splice annotation | Index files; preparation is separate from recording |
| `capzyme index-record` | Existing FASTA/index | Dictionary and checksum manifest; no index building |
| `capzyme trim` | Raw single-end FASTQ | Cutadapt-trimmed FASTQ and QC |
| `capzyme rna-align` | Trimmed FASTQ and RNA FASTA/index | Bowtie2 alignments, accepted BAM, genuinely unmapped FASTQ and accounting |
| `capzyme genome-align` | Trimmed or RNA-unmapped FASTQ and genome FASTA/index | Bowtie2/HISAT2 alignments, decoy/CIGAR-filtered accepted BAM and accounting |
| `capzyme count` | Ordered single-end BAMs and matching reference | Original integer, one-based single-nucleotide matrix using the shared BEDTools engine |
| `capzyme bigwig` | Unthresholded matrix or BAMs | Separate strand tracks and factors; shared counting and normalization engines |

`workflow/scripts/compare_matrices.py` provides exact keyed regression comparison;
`validate_bigwigs.py` independently checks track values. These are validation
utilities rather than analysis stages. See [coordinates and counting](tools/counting.md),
[genome policy](genome_mapping.md) and [BigWig methods](tools/bigwigs.md).

## Snakemake orchestration

| Entry point | Coordinated stages |
| --- | --- |
| `workflow/Snakefile` | RNA FASTQ trimming → generated/fixed reference → Bowtie2 → counts; optionally sequential genome mapping/counts |
| `workflow/genome.smk` | Already trimmed FASTQ → existing genome index → genome counts; no trimming |
| `workflow/existing_bam.smk` | Read-only existing BAMs → output-local filtering/ordering as needed → same counter; no mapping |
| Optional BigWig targets | Unthresholded stage matrices → standalone BigWig engine; RNA and genome normalization stay separate |

The launcher (`capzyme workflow`) validates configuration and schedules the DAG.
Rules expand samples, manage dependencies/resources and publish provenance.
They do not implement separate counting or normalization algorithms. Version-1
workflow configuration and the optional version-2 standalone settings adapter
remain distinct supported interfaces; see [settings](tool_settings.md).

`core.py` and `bam_count.py` implement shared BAM filtering and BEDTools
`genomecov -5 -dz` assembly, with one-based matrix coordinates. `run_step.py`,
`bam_step.py`, `genome_step.py` and `bigwig_step.py` are orchestration adapters.
`workflow_bigwig.py` validates cohort settings and publishes outputs safely.
The [reference contract](reference_verification.md) distinguishes automatic
observed hashes from optional strict trusted expectations.

BigWig raw/RPM/supplied modes require no R. Optional DESeq2/control modes invoke
actual R/DESeq2 functions in `R/deseq2_normalization.R`. Future independent
statistical analysis will consume the original integer matrix, never track values;
neither analysis nor BigWig execution is a prerequisite of the other.
MultiQC orchestration and reusable differential-expression analysis remain deferred.

See the [complete workflow examples](../examples/workflows/README.md) and
[current finalization evidence](finalization_report.md). Historical experiment
reports are retained separately for owner/publication review, not installation instructions.
