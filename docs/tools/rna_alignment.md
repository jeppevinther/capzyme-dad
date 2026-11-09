# Align reads to an RNA reference with Bowtie2

**Program:** Bowtie2 RNA-reference alignment with CapZyme filtering, accounting and
handoff verification. Executable adapter:
[rna_handoff.py](../../workflow/scripts/rna_handoff.py). Requires Python, Bowtie2
mapper/inspector, SAMtools, trimmed single-end FASTQ, unchanged finished RNA FASTA
and complete existing Bowtie2 index. Collapse is never implicit.

```bash
capzyme rna-align --input-file "$TRIMMED_FASTQ" \
  --reference "$RNA_FASTA" --index "$RNA_INDEX" --seed 0 --mapq 20 \
  --threads 4 --output-dir "$WORK/new-rna-alignment"
# Optional strict settings/manifest route:
capzyme rna-align --settings "$WORK/tools.yaml" \
  --sample A --output-dir "$WORK/new-rna-alignment" --preview
```

Policies: `--very-sensitive`, default end-to-end; MAPQ≥20, excluded flags2308,
retain duplicate/QC-fail records; seed0. `--norc` forbids reverse-complement
alignments relative to each reference record. The direct script requires this
flag explicitly; the installed command and settings adapter default
`allow_reverse_complement: false` and prints `--norc`. Either supports enabling
reverse complements deliberately (omit direct `--norc`, or override the adapter
with `--allow-reverse-complement`). No alternative preset is advertised.

Native input production is also supported:

```bash
bowtie2 --very-sensitive --end-to-end --norc --seed 0 -p 4 \
  -x "$RNA_INDEX" -U "$TRIMMED_FASTQ" --un "$UNMAPPED_FASTQ" -S "$ALL_SAM"
samtools view -b -q 20 -F 2308 "$ALL_SAM" | samtools sort -m 256M -o "$ACCEPTED_BAM"
```

The adapter's exact executed native commands and working-directory aliases are
in `commands.jsonl`. It additionally verifies dictionaries, source hashes,
unique primary read identity and unmapped sequence/quality completeness.
Only aligner-unmapped reads enter `unmapped.fastq`; low-MAPQ mapped rejects do
not. Counts are alignment records, not deduplicated molecules or paired fragments.

Exact outputs: `accepted.bam` (filtered, sorted), `alignments.bam` (all records,
sorted), `unmapped.fastq`, `qc.json`, `provenance.json`, `commands.jsonl`,
`reference.sizes.tsv`, `fasta.sizes.tsv`, `index.summary.txt`; toolbox additionally
writes `methods.json`. No BAM index is required for standalone counting; create
a new managed index only if a browser needs one. Use a new external output
directory. Native BAMs can go directly to the counter, with the reduced
accounting guarantees described in [the toolbox](../toolbox.md).
