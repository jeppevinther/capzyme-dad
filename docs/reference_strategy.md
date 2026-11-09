# Generated RNA reference strategy

The implemented reference space is a collapsed ncRNA FASTA plus zero or more
explicitly supplied viral FASTAs. No reference is downloaded automatically.
Record source release, assembly, accessions, orientation and checksums in the
experiment's private reference records. Run provenance hashes the actual inputs.

## Established collapse semantics

The workflow imports the unchanged, checksum-verified `merge_by_prefix_v2.py`
and passes `k=30` explicitly. Its standalone default remains 70. A configuration
can choose k=70 for a separate comparison run without modifying the original.

- Whitespace is removed and sequences are uppercased before grouping.
- The group key is the first k bases; sequences shorter than k use their entire
  sequence. Ambiguous characters in keys are matched literally, not as wildcards.
- The representative is the longest full-length sequence in the group.
  Equal-length representatives retain the first input record. Input order is
  therefore part of provenance; reordering can change representatives and IDs.
- Group order follows first encounter. The most common family/symbol labels the
  group, with first-encounter tie handling inherited from Python's Counter.
  Original family simplification, safe-name truncation, member counts and
  collision suffixes are retained. No sequence-based tie breaker is introduced.
- Description fields mentioning the word pseudogene are always removed by the
  original algorithm. Configurable named biotypes are also removed; no named
  biotype exclusion is silently imposed on every experiment. Empty sequences
  are skipped as in the original. Run logs retain record and filter summaries.
- `reference/ncrna.tsv` preserves the original membership-table columns,
  representative ID/length, group prefix, member IDs and symbols. No membership
  reassignment or genomic lift-over is inferred.

The assembled reference is validated for nonempty IUPAC DNA sequences and unique
record IDs. An invalid/empty combined reference fails before indexing. Viral
records are preserved in provided order before collapsed ncRNAs, matching the
original experiment's assembly order. They are not collapsed. Explicitly supplied
plus- and minus-sense viral records remain separate reference entries.

## Why k=30, and what comparison with k=70 means

Shared 5′ prefixes in multicopy short RNAs can create numerous nearly equivalent
reference entries. Collapsing them reduces reference complexity and can reduce
the dispersion of accepted starts across copies, increasing counts and power for
the resulting group-level feature. A count assigned to that group is not evidence
for a particular original gene copy.

A shorter prefix merges more references. This may reduce multimapping among
members, but it also removes alternative downstream sequences. Retaining only
the longest full-length sequence can lose mappings of reads that match another
member beyond the shared prefix, or change MAPQ through changes in competing
alignments. Thus k=30 cannot be assumed to improve sensitivity or power for every
RNA. k=70 can preserve more distinguishable sequences but retain more mapping
ambiguity, particularly for reads shorter than the distinguishing region.

An eventual k=30/k=70 comparison should hold input references, read processing,
aligner/version, filtering and sample set constant. Compare group/member counts,
retained reads, MAPQ distributions, ambiguous mappings, observed 5′ sites, and
downstream feature-level statistical behavior. Use separate output roots and
reference IDs/membership tables; do not directly equate rows with different
group definitions. No real-data comparison or alternative collapse algorithm is
implemented by this reference-preparation interface.

## Orientation and interpretation

With `--norc`, Bowtie2 aligns only the input-read orientation to each supplied
reference. A read matching a separately supplied minus-sense viral reference
can therefore have SAM alignment strand `+` while representing viral minus-sense
RNA. SAM strand is relative to that reference sequence, not an automatic label
of biological viral polarity. Keep reference orientation metadata explicit.

The matrix's chromosome field is a reference FASTA ID, and its positions refer
to that entry. A 5′ position in a collapsed ncRNA is a transcript/group coordinate,
not a genomic coordinate. On a full-length reverse-complement viral record of
length L, a one-based base p corresponds to L-p+1 on its opposite sequence;
that simple transform applies only when records are known full reverse
complements. No such conversion is performed implicitly.

The governing mapping behavior is documented in the official
[Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
Reference grouping itself is governed by the supplied script, not by a new
algorithm inferred from documentation.
