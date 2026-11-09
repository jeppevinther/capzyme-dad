# Collapse an RNA reference

**Program:** unchanged [merge_by_prefix_v2.py](../../legacy_scripts/merge_by_prefix_v2.py).
The optional safe CLI below calls `core.collapse`, the same hash-verified function
used by Snakemake. Python standard library only; no Snakemake dependency.

```bash
capzyme collapse --reference "$NCRNA_FASTA" \
  --prefix-length 30 --output-dir "$WORK/new-collapsed-reference"
# Optional, repeatable: --drop-biotype lncRNA
```

Input is an ncRNA FASTA with the intended reference-relative orientations.
The original algorithm groups equal 5′ prefixes, chooses the longest sequence,
breaks equal-length ties by input order, and preserves membership labels. Prefix
length and biotype exclusions are scientific choices, not mapping parameters.
Review [reference strategy](../reference_strategy.md) for short/pseudogene records.

Outputs: `collapsed.fa` (representatives), `collapsed.tsv` (membership),
`methods.json` (settings, reference/legacy hashes and grouping diagnostics).
The output directory must be new/external; source stays unchanged. CapZyme adds
source verification, safe output creation and provenance, not a second collapse
algorithm. Combining viral and representative records is a separate explicit
reference preparation step; general mapping never collapses its input FASTA.
