# Trim FASTQ with Cutadapt

**Program:** Cutadapt; optional adapter `capzyme trim`. Supply raw, single-end
Phred33 FASTQ, adapter and reviewed trimming settings. Python plus Cutadapt;
Snakemake is unnecessary.

```bash
capzyme trim --input-file "$RAW_FASTQ" \
  --adapter "$ADAPTER" --quality 20 --minimum-length 15 --nextseq-trim 20 \
  --threads 4 --output-dir "$WORK/new-trimming"
```

Defaults are quality 20, minimum length 15, NextSeq threshold 20 and one thread.
The adapter is required; these defaults are documented behavior, not universal
library recommendations. It executes:

```bash
cutadapt -a "$ADAPTER" -q 20 -m 15 --nextseq-trim 20 --cores 4 \
  --json "$WORK/new-trimming/cutadapt.json" \
  -o "$WORK/new-trimming/trimmed.fastq.gz" "$RAW_FASTQ"
```

Outputs: trimmed FASTQ, Cutadapt's JSON/log and `methods.json` recording resolved
command, version and status. CapZyme adds settings resolution and safe new output
selection. Cutadapt performs all trimming. Use the retained FASTQ for mapping;
do not label original raw reads as already trimmed. Native Cutadapt outputs are
also supported inputs to the downstream tools.

## Disable extra two-color trimming

Use `--no-nextseq-trim`, or `trimming.nextseq_trim: null` in standalone settings.
In workflow settings use `preprocessing.nextseq_trim: null`. This omits the
Cutadapt option entirely and keeps ordinary quality trimming active. Numeric
zero still enables the additional treatment and can trim high-quality trailing
Gs in Cutadapt 4.9. Nonzero behavior is unchanged. CLI numeric/off choices are
mutually exclusive and override settings, including explicit null; omitted
standalone options/settings default to 20. Workflow preprocessing fields remain
explicit. Choose this treatment for suitable two-color NextSeq/NovaSeq data;
other chemistries can use ordinary adapter/quality trimming.
