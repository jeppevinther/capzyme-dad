# Align reads to a genome with HISAT2 or Bowtie2

**Programs:** HISAT2 for splice-aware mapping; Bowtie2 for unspliced mapping.
Adapters are [hisat2_genome.py](../../workflow/scripts/hisat2_genome.py) and
[bowtie2_genome.py](../../workflow/scripts/bowtie2_genome.py), sharing acceptance,
accounting and provenance functions. They are not new alignment algorithms.
Requires the selected mapper/inspector, SAMtools, Python, trimmed single-end FASTQ,
matching FASTA/existing index and explicit decoy policy.

```bash
capzyme genome-align --aligner hisat2 \
  --input-file "$TRIMMED_OR_RNA_UNMAPPED" --reference "$GENOME_FASTA" \
  --index "$HISAT2_INDEX" --decoys "$DECOY_IDS" --mapq 20 --seed 0 \
  --threads 4 --output-dir "$WORK/new-genome-alignment"
# Or select --aligner bowtie2 with a Bowtie2 index for an unspliced genome.
# For a reference with no decoys, explicitly substitute --no-decoys.
capzyme genome-align --settings "$WORK/tools.yaml" \
  --sample A --output-dir "$WORK/new-genome-alignment" --preview
```

HISAT2 enforces `--very-sensitive --no-softclip`; Bowtie2 enforces
`--very-sensitive --end-to-end`. Both orientations are searched, seed0, MAPQ≥20,
exclude flags2308, retain duplicates/QC-fail by default. Explicit decoys are
excluded. Conservative CIGAR acceptance requires terminal M/=/X, no S/H/P;
internal I/D/N are allowed. `--uncertain-cigar reject|error` selects handling of
uncertain endpoints. Duplicate exclusion optionally adds flag1024. Settings
reject alternative presets/soft-clipping modes. Native examples:

```bash
hisat2 --very-sensitive --no-softclip --seed 0 -p 4 -x "$HISAT2_INDEX" \
  -U "$TRIMMED_FASTQ" --un "$UNMAPPED_FASTQ" -S "$ALL_SAM"
bowtie2 --very-sensitive --end-to-end --seed 0 -p 4 -x "$BOWTIE2_INDEX" \
  -U "$TRIMMED_FASTQ" --un "$UNMAPPED_FASTQ" -S "$ALL_SAM"
```

Native alignment alone does not apply the CapZyme decoy/CIGAR policy. Generic
SAMtools filtering plus standalone counting cannot infer missing evidence or
retroactively supply the adapters' full accounting guarantees.

Outputs: `accepted.bam`, `alignments.bam`, `unmapped.fastq`, `qc.json`,
`provenance.json`, `commands.jsonl`, `mapping_by_reference.tsv`,
`reference.sizes.tsv`, `fasta.sizes.tsv` (when FASTA supplied), `index.summary.txt`;
toolbox adds `methods.json`. Input identity/qualities and exact unmapped handoff
are checked; source inputs/indexes are rehashed unchanged. Sort memory defaults
to 256MiB per worker. New settings always require an existing index; the legacy
Bowtie2 direct CLI's optional map-and-build behavior is preserved but not the
recommended settings route. [Detailed guide](../genome_mapping.md) and
[separate preparation](index_preparation.md).
