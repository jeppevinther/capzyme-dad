# Complete editable workflow examples

Copy one sample sheet and its settings file into your external analysis folder:

- [Sequential samples](sequential.samples.tsv) and [sequential settings](sequential.settings.yaml): raw FASTQs → Cutadapt → fixed RNA Bowtie2 mapping → genuinely aligner-unmapped reads → HISAT2 genome mapping.
- [Direct-genome samples](genome.samples.tsv) and [direct-genome settings](genome.settings.yaml): already trimmed FASTQs → Bowtie2 genome mapping. Trimming is not part of this route.

Replace placeholder paths, provide your adapter for the sequential route, and
review chemistry, decoys, threads and memory. `nextseq_trim: null` uses ordinary
quality trimming; choose a numeric threshold for suitable two-color data.
Zero is still an enabled numeric setting. Do not forward mapped RNA rejects to
the genome as though they were aligner-unmapped.

RNA-first assignment is useful when prepared RNA records should receive a
priority assignment before competing genomic copies. It is a scientific choice:
reference composition changes assignment, and direct genome mapping is a valid
alternative for appropriate experiments, including some eukaryotic analyses.

Expected hashes are optional; observed hashes and compatibility checks are
automatic. Retain trustworthy native build records for historical identity and
annotation provenance. Both examples generate raw/RPM tracks from unthresholded
counts; the workflow computes complete accepted depths separately for RNA and
genome, with one sample factor shared across genomic strands. Optional DESeq2
blocks require a separate R environment and an explicitly selected estimator.

Run `capzyme workflow check`, `preview`, then `run` with `--route sequential`
or `--route genome`, following the [README](../../README.md#snakemake-workflow).
These are templates, not bundled research data or ready-to-run biological analyses.
