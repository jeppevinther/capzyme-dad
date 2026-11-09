# Standalone genome mapping

Two independent programs map single-end **already trimmed** Phred33 FASTQ:
`workflow/scripts/bowtie2_genome.py` for unspliced genomes and
`workflow/scripts/hisat2_genome.py` for splice-containing genomes. Neither runs
ncRNA mapping, trimming, counting or normalization internally. Existing counting
and BigWig algorithms are unchanged. Both accept `--threads`, new external
`--output-dir`, `--reference-sha256`, `--mapq`, `--sort-memory`,
`--exclude-duplicates`, `--seed` and `--uncertain-cigar reject|error` where relevant.
Run either executable with `--help` and `map --help` / `prepare --help` for inputs.

## Bowtie2 Genome

```bash
python workflow/scripts/bowtie2_genome.py map \
  --fastq "$TRIMMED_FASTQ" --reference "$GENOME_FASTA" \
  --threads 4 --output-dir "$NEW_OUTPUT"

# Explicit reusable preparation; genome-only analysis needs no ncRNA stage.
python workflow/scripts/bowtie2_genome.py prepare \
  --reference "$GENOME_FASTA" --threads 4 --output-dir "$NEW_INDEX_OUTPUT"
python workflow/scripts/bowtie2_genome.py map \
  --fastq "$TRIMMED_FASTQ" --index "$NEW_INDEX_OUTPUT/index/genome" \
  --reference "$GENOME_FASTA" --threads 4 --output-dir "$NEW_OUTPUT"
```

The fixed alignment policy is `--very-sensitive --end-to-end`, with both genomic
orientations enabled. No GFF/GTF is needed. Multiple chromosomes, plasmids and
other sequences remain separate dictionary entries. Existing complete `.bt2` or
`.bt2l` families are supported; mixed/incomplete families are rejected. Supplying
FASTA alone builds an output-local index. See the
[Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

## HISAT2 Genome

```bash
python workflow/scripts/hisat2_genome.py map \
  --fastq "$TRIMMED_FASTQ" --index "$HISAT2_PREFIX" \
  --reference "$GENOME_FASTA" --decoys "$DECOY_IDS" \
  --threads 4 --output-dir "$NEW_OUTPUT"

# Separate explicit preparation, illustrated for small validation references.
python workflow/scripts/hisat2_genome.py prepare \
  --reference "$GENOME_FASTA" --gtf "$COMPATIBLE_GTF" \
  --decoys "$DECOY_IDS" --threads 4 --output-dir "$NEW_INDEX_OUTPUT"
```

The fixed policy is `--very-sensitive --no-softclip`, both orientations enabled;
there is no RNA-library strandedness restriction. Mapping requires an existing
complete `.ht2` or `.ht2l` family. Preparation optionally accepts `--gtf`, or
pre-extracted `--splice-sites` / `--exons` tables (not both routes together).
GTF exon coordinates are 1-based inclusive; extracted tables use the HISAT2
four-column zero-based inclusive format: reference, left, right, strand. Unknown
identifiers, out-of-bounds positions, absent exon records or malformed tables
fail. The installed graph builder fails on exon-only annotations without splice
sites. This tool rejects that case explicitly (including a single-exon-only GTF);
it never silently discards annotations. For a linear index, explicitly prepare
without annotations. Splice sites alone or paired splice/exon tables are supported. Names and
coordinate bounds are checked; this does not prove annotation/assembly biological
compatibility. The caller must select the same assembly/release.

Either `--decoys` or `--no-decoys` is required for HISAT2 prepare **and** map.
The list is headerless, one exact unique reference identifier per line, with no
heuristic guessing. Reuse the same designation when mapping a prepared index;
provenance records both decisions. A genuine no-decoy reference can use
`--no-decoys`. For historical GRCh38 `no_alt_plus_hs38d1`, keep all reference and
decoy sequences in the original index/dictionary and provide its verified decoy
list. No actual human index or decoy inventory has been audited in this milestone.

Annotation-aware human index construction can require roughly 200 GB RAM per the
[HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/); the small synthetic
fixture does not benchmark that cost. Four threads or sort-memory settings do
not cap index-builder memory. Full human index preparation is outside this work.

The shared code uses the matching native HISAT2 inspector (`-s`/`-l`) when
available, avoiding an old vendor wrapper's Python-`imp` incompatibility. HISAT2
mapping uses read-only output-local symlinks and relative arguments so its wrapper
also works with whitespace in source/output paths. The wrapper still implements
`--un`; it is never patched. Effective commands, alias targets and working
directory are recorded. No installed package was modified for these workarounds.

## Output and coordinate/filter policy

| Output | Meaning |
| --- | --- |
| `alignments.bam` | All original aligner records, sorted; includes unmapped, secondary, decoy and filter-rejected records |
| `accepted.bam` | Sorted primary, non-decoy, coordinate-certain records passing the declared filtering policy; use this for genomic counting |
| `unmapped.fastq` | Exact original sequence/quality for reads the aligner classified unmapped; order is not an API guarantee |
| `reference.sizes.tsv` | Index reference identifiers and lengths, including decoys |
| `qc.json` | Exclusive primary-read outcomes, total/nonprimary records, decoy and handoff counts |
| `mapping_by_reference.tsv` | All dictionary entries with lengths, decoy labels and mapped primary counts |
| `commands.jsonl` | Command/status/timing JSON interleaved with raw tool logs; not a pure JSONL stream |
| `provenance.json` | Completion/failure, source/index hashes, policies, code/executable identities, versions and output hashes |

All primary input read names must be unique (first whitespace-delimited token).
Strict four-line FASTQ and gzip are supported; multiline FASTQ/paired reads are
not. Disk-backed bookkeeping checks each input has exactly one primary SAM
record and that forwarded sequences/qualities are exactly the original unmapped
subset. This prevents mapped-but-rejected reads from being forwarded. SQLite
bookkeeping and temporary SAM/sort files require output-local disk space and have
not been benchmarked on full genomes. Inputs are hashed before/after; no source
FAI, BAM index or other source-side file is written. Existing outputs are refused.

Default filtering reuses the existing policy validator: MAPQ >=20; exclude flags
2308 (unmapped, secondary, supplementary); retain duplicates and QC-failure
records. `--exclude-duplicates` additionally excludes flag 1024. In the accepted
BAM, declared decoys are also excluded. Exclusive mapped rejection reasons are
flags, then decoy, then uncertain CIGAR, then low MAPQ. Per-reference mapped
primary statistics precede these acceptance filters. Rejected mapped reads stay
in `alignments.bam`; they are never automatically remapped. Secondary/supplementary
records are separately accounted for and never counted as extra input reads.

For coordinate certainty, terminal operations must be M, = or X. Any S/H/P and
terminal I/D/N are rejected by default, even if only the other end is affected;
this deliberately conservative rule is more restrictive than the generic counter.
`--uncertain-cigar error` instead fails on such primary mapped records. Internal
I/D and splice N are permitted, with reference bounds and query consumption
validated. No alignment is repaired or coordinate-shifted. The unchanged BEDTools
method then counts plus at the leftmost and minus at the rightmost aligned
nucleotide, including intron spans for minus-strand endpoints. Matrices retain
one-based `start=end`; BigWigs retain zero-based `[p-1,p)`.

## Counting, sequential use and later Snakemake orchestration

```bash
python workflow/scripts/count_5ends.py "$GENOME_OUT/accepted.bam" \
  --reference "$GENOME_FASTA" --sample-name sample --strand both \
  --output "$NEW_COUNTS_TSV"
python workflow/scripts/bigwig2.py --bam "$GENOME_OUT/accepted.bam" \
  --reference "$GENOME_FASTA" --sample-name sample --cohort sample \
  --normalization raw --strands both --output-dir "$NEW_TRACKS_DIR"
```

Always pass `accepted.bam`, not the all-records archive, to preserve the genome
CIGAR/decoy policy. The counter re-applies its established MAPQ/flag policy without
modifying inputs. In index-only mapping, supply the matching original FASTA to
the downstream tools; no reconstructed genome FASTA is implicitly generated.
If both index and FASTA are given, identifiers/lengths are checked, but dictionary
agreement alone does not prove nucleotide identity. All individual input hashes
are recorded; reference selection remains the caller's explicit responsibility.

For a human/viral experiment, the first ncRNA/viral aligner must export its own
aligner-unmapped FASTQ. Pass exactly that file to HISAT2 Genome. If the first
stage has no such artifact, stop and review its extraction rather than feeding
all MAPQ-rejected reads. Synthetic validation exercises this handoff using the
standalone Bowtie2 tool against a tiny artificial ncRNA space, then HISAT2; it
does not claim the established Snakemake RNA stage already exposes this output.
Maintain separate ncRNA and genome matrices and separate reference dictionaries.

Future Snakemake rules should declare complete index families, FASTQ, optional
FASTA/annotations and decoy lists as inputs, call these independent programs,
and expose accepted BAM/QC/provenance/unmapped artifacts. Use a fresh per-attempt
directory and publish only a successful complete output; rules must account for
the CLIs' refusal to reuse existing directories. Subsequent rules call the
unchanged counter and BigWig CLI. CPU and memory scheduling belong to orchestration
and external execution policy; this milestone adds no active workflow rule.

Example [human](../examples/genome_mapping/human.yaml) and
[bacterial/archaeal](../examples/genome_mapping/bacterial_archaeal.yaml) YAMLs
are CLI argument bindings for future orchestration, not a new `--config` loader.
Shared normalization across reference spaces remains unresolved: equal supplied
factors require matching sample identities, cohort and compatible quantitative
scope. Genomic-only RPM denominators do not equal ncRNA-stage denominators, and
saved BigWig factor provenance may reject another dictionary. Do not bypass those
checks or combine normalized values as raw counts. No new normalization method,
internal-background approach or cross-space factor transfer was implemented.

## Reference preparation and compatibility

Inspector descriptions are separated from first-token sequence IDs to match
FASTA/SAM dictionaries. A dictionary match alone cannot establish the sequence
or embedded annotation origin; retain build records and use optional trusted
manifests as described in [reference verification](reference_verification.md).

For annotated HISAT2 preparation, reconcile GTF sequence identifiers with the
exact FASTA dictionary, verify assembly/strand/coordinate compatibility and
explicitly account for absent sequences. Do not guess accession aliases. Native
`hisat2_extract_splice_sites.py` and `hisat2_extract_exons.py` generate tables for
`hisat2-build --ss ... --exon ...`. Check nonzero tables, valid intervals and
embedded annotation counts after building, in addition to all eight index parts
and the reference dictionary. Build into a new managed directory with a reviewed
human-scale resource budget; never overwrite a historical index.

HISAT2 can also accept known splice sites at mapping time, but the current
standalone mapping CLI does not expose that option. Prefer a prepared annotated
index for routine splice-aware mapping. See [native preparation commands](tools/index_preparation.md)
and [optional sequential workflow](sequential_mapping.md).
