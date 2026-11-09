# Count single-nucleotide 5′ ends

**Programs:** SAMtools filtering/order checks and BEDTools `genomecov -5 -dz`.
[count_5ends.py](../../workflow/scripts/count_5ends.py) and optional toolbox settings
both call the validated `bam_count` / `core.matrix` engine used by Snakemake.
No alternative counting algorithm is introduced.

```bash
capzyme count --reference "$FASTA" \
  --bam "$BAM_A" --sample-name A --bam "$BAM_B" --sample-name B \
  --strand both --mapq 20 --output-dir "$WORK/new-counts"
capzyme count --settings "$WORK/bam-tools.yaml" --preview
```

Required: one or multiple single-end BAMs and a matching unchanged FASTA; explicit
ordered sample names are recommended. Default: MAPQ≥20, flags2308 excluded,
duplicates and QC-failure reads retained, both strands, no total-count threshold.
`--exclude-duplicates`, `--strand plus|minus|both` and `--min-total-count N` expose
the supported choices. Each eligible alignment contributes one aligned 5′ end.
Plus uses the left aligned base; minus uses the right reference-consuming end.
This generic counter accepts CIGARs already present; apply the genome adapter's
conservative acceptance first when that scientific policy is required.

Sources are never modified/reindexed/resorted. Filtering uses temporary copies;
sorting is skipped if those copies are actually ordered. Already discarded reads
cannot be recovered. Output TSV: `chromosome`, `strand`, 1-based `start`, `end`
(equal), then ordered integer sample columns. Direct CLI outputs exactly the
requested matrix and default `OUTPUT.jsonl` provenance/command log. A nonzero
threshold filters this requested standalone matrix; **retain a separate
unthresholded matrix for BigWigs/statistics**. The optional workflow preserves
raw plus an extra thresholded matrix.

The toolbox emits `counts.tsv`, `counts.tsv.jsonl`, `methods.json` in a new
external directory. Details and BAM-input limitations:
[existing-BAM counting](../existing_bam_counting.md).
