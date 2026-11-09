Start with the [tool catalogue](../../docs/toolbox.md). `tool_settings.py` backs
the optional standalone `toolbox.py` adapter, reusing these scientific engines.
`resource_guard.py` is Linux operational supervision with independent storage
scans; `validate_bigwigs.py` independently checks decoded tracks against integer
matrices using disk-backed SQLite. Neither is a second scientific counter.

Milestone 1 helpers:

- configuration.py validates the implemented YAML/manifest contract.
- core.py wraps unchanged reference collapse, validates/combines FASTAs, assembles
  BEDTools 5′ counts and applies optional matrix row thresholds.
- run_step.py connects Python helpers to Snakemake and writes run provenance.

Original reference scripts remain unchanged in `legacy_scripts/`. The original
collapse module is imported only after checksum verification; its grouping and
representative-selection algorithm are not reimplemented.

Milestone 2 preparation adds fixed_reference.py: checksum and ordered-content
validation for immutable external FASTA/index inputs. Configuration calls its
hash checks during DAG preparation; the fixed-reference rule repeats them and
inspects index content before allowing any FASTQ preprocessing.

`pilot_validation.py` implements the completed, explicitly authorized bounded
old/new comparison. It takes private runtime paths, never chooses a host, fixes
the four-library/100,000-record boundary, inherits nice=10, restricts affinity,
monitors host/resources and retains outputs only in a fresh server run directory.
It is not automatically called by the scientific workflow. Review and obtain
new approval before rerunning it or expanding its scope.


Milestone 3A adds independently executable `count_5ends.py` and
`compare_matrices.py` (`python SCRIPT --help`). `bam_count.py` shares BAM
preparation and the existing matrix/threshold functions; `bam_configuration.py`
and `bam_step.py` connect these tools to `workflow/existing_bam.smk`. No second
scientific counting algorithm is introduced. Actual-order validation skips
redundant sorting and all indexing. See `docs/existing_bam_counting.md`.


Milestone 4A adds standalone `bigwig2.py` and `bigwig_normalization.py`, with an
internal `deseq2_factors.R` adapter. Actual DESeq2 estimation lives in the reusable
`R/deseq2_normalization.R` function. Raw/RPM/supplied modes do not invoke R.
No new Snakemake target is activated. See `docs/bigwig2.md` and the validation
report for all five modes and synthetic coverage.

Milestone 5A adds two separate genome CLIs: `bowtie2_genome.py` for unspliced
genomes and `hisat2_genome.py` for splice-aware genomes. `genome_common.py` shares
input/index validation, BAM/QC handling and provenance. Both produce sorted
accepted BAMs and exact aligner-unmapped handoffs; counting/normalization remain
in their unchanged standalone programs. See `docs/genome_mapping.md`.

Sequential orchestration uses `rna_handoff.py` for the established RNA Bowtie2
policy plus exact unmapped-FASTQ preservation. `genome_step.py` only adapts
Snakemake inputs/resources to the standalone CLIs and shared counter;
`sequential.py` publishes fresh attempts and verifies exclusive read identities.
`genome_configuration.py` validates external references, complete index families
and policy dependencies. See `docs/sequential_mapping.md` for inputs, outputs,
restart semantics and genome-only operation.
