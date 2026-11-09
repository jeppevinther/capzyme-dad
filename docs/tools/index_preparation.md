# Build and verify a reference index

**Programs:** native Bowtie2/HISAT2 builders and inspectors. No wrapper is required
to build an index. Use a finished FASTA, appropriate aligner, CPU/storage budget
and a fresh external directory; never overwrite a historical prefix.

```bash
INDEX_DIR="$(mktemp -d "$WORK/bowtie2-index.XXXXXX")"
bowtie2-build --threads 4 "$FASTA" "$INDEX_DIR/reference"
bowtie2-inspect -s "$INDEX_DIR/reference" > "$INDEX_DIR/summary.txt"
capzyme index-record --aligner bowtie2 \
  --reference "$FASTA" --index "$INDEX_DIR/reference" \
  --output-dir "$WORK/new-bowtie2-verification"
```

Bowtie2 emits six `.bt2` or `.bt2l` parts. For an appropriately reviewed splice-aware
build, use compatible annotation (same assembly and exact sequence identifiers):

```bash
INDEX_DIR="$(mktemp -d "$WORK/hisat2-index.XXXXXX")"
hisat2_extract_splice_sites.py "$COMPATIBLE_GTF" > "$INDEX_DIR/splice_sites.tsv"
hisat2_extract_exons.py "$COMPATIBLE_GTF" > "$INDEX_DIR/exons.tsv"
hisat2-build -p 4 --ss "$INDEX_DIR/splice_sites.tsv" \
  --exon "$INDEX_DIR/exons.tsv" "$FASTA" "$INDEX_DIR/reference"
capzyme index-record --aligner hisat2 \
  --reference "$FASTA" --index "$INDEX_DIR/reference" \
  --output-dir "$WORK/new-hisat2-verification"
```

HISAT2 emits eight `.ht2`/`.ht2l` parts. This example is not a human-build resource
request: splice-annotated human builds need separate memory/storage/runtime review.
Omitting annotation deliberately builds a linear index. Existing managed
`bowtie2_genome.py prepare` and `hisat2_genome.py prepare` remain optional; the
[full guide](../genome_mapping.md) covers their GTF/table checks and outputs.

`index-record` reuses existing complete-family and inspector/dictionary checks.
It writes `reference_manifest.json` (FASTA hash, all suffix hashes, dictionary,
aligner and trust statement), `reference.sizes.tsv`, `fasta.sizes.tsv`,
`index.summary.txt`, `commands.jsonl` and `methods.json`. Native HISAT2 inspectors
are selected where the vendor Python wrapper is incompatible with Python 3.12.
Original FASTA/indexes remain read-only; mixed or incomplete families fail.

**Trust distinction:** first registration records newly observed hashes and
identifier/length agreement. It does not prove historical sequence identity or
embedded annotation correctness. Retain trustworthy builder/reference provenance;
use `--trusted-manifest "$PREVIOUS_MANIFEST"` to require existing pins before
recording. Never replace a mismatching trusted hash with a new observation.
Relocatable suffix manifests allow reuse without hand-entering every part hash.
The version-1 fixed RNA workflow additionally reconstructs index sequence content;
the new generic registration does not claim that stronger check.

Settings mapping accepts only `action: use_existing` and verifies its manifest.
`build_new` is rejected with preparation guidance. It never repairs a missing
index automatically. Reference annotations/decoys and any future changed index
belong in distinct managed locations.
