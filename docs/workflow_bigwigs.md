# Optional workflow BigWigs

The RNA, sequential RNA-to-genome, genome-only and existing-BAM workflows accept
an optional top-level `bigwigs` mapping in their **version-1** configuration.
It is separate from the standalone version-2 `bigwig` schema. Omission schedules
no tracks and adds no pyBigWig/R dependency. See [example](../config/bigwigs.example.yaml).

Each reference space (`rna`, `genome`) has `methods`, optional `cohort`,
`emit_samples`, `strands` and `mem_mb` (default 8192). Only applicable spaces are
allowed. Cohort defaults to all configured sample IDs in sample-sheet order;
emitted samples default to that cohort. Explicit cohorts and emitted subsets
must contain unique IDs; emitted samples must belong to the cohort. Selected
strands must exist in the matrix: a plus-only matrix cannot produce minus tracks.

| Method | Options inside the method mapping |
| --- | --- |
| `raw` | `{}`; original integer signal |
| `rpm` | `{}` audits complete accepted BAM record totals; optional `depths` replaces that audit |
| `deseq2` | Required `deseq2_type: ratio` or `poscounts`; optional `rscript` |
| `controls` | Required `deseq2_type`, `controls` TSV; optional `min_control_features` (10), `rscript` |
| `supplied` | Required `size_factors` TSV and `factor_kind: size-factor` or `multiplier`; optional `factor_metadata` |

Path fields follow workflow semantics: relative to the external working
directory (configuration directory with the launcher). All control/factor files
must be outside the output directory. DESeq2/control modes probe the selected
Rscript for DESeq2 before scheduling jobs. Unknown fields, irrelevant method
options, missing files, invalid factor tables, malformed controls and invalid
cohort selections fail at configuration time. Whether controls are present and
informative in newly produced counts is tested by the unchanged scientific
engine after counting; failures do not trigger an estimator switch.

RNA factors and genome factors are always estimated separately. A genome sample
has one factor shared by its plus and minus tracks. The input is always
`counts/5ends.raw.tsv` or `counts/genome.5ends.raw.tsv`, never a thresholded table.
The workflow does not recount 5′ ends for tracks. For automatic RPM it runs
`samtools view -c` on accepted BAMs and verifies each BAM's checksum before/after,
saving `bigwigs/SPACE/accepted_depths.tsv` and `logs/bigwig.SPACE.depths.log`.
The denominator includes all accepted references and both alignment strands,
even if the count matrix or displayed signal is plus-only. External depths need
`sample_id`, positive `accepted_records`, and one nonempty `scope_id` describing
that same complete accepted pool; the caller is responsible for the audit origin.
Empty libraries cannot support RPM or DESeq2; choose a scientifically appropriate
explicit cohort or use raw tracks. No silent removal of empty samples occurs.

## Outputs and recovery

Each `bigwigs/SPACE/METHOD/` contains:

- `tracks/SAMPLE.plus.bw` and/or `tracks/SAMPLE.minus.bw`.
- `normalization_factors.tsv`, `reference.sizes.tsv`, `provenance.json`.
- `commands.log` (native engine calls) and `workflow.json` (exact Python command,
  attempt location, exit status).
- R modes: `normalization/estimation_counts.tsv`, `r_size_factors.tsv`,
  `r_diagnostics.tsv`, `r_diagnostics.tsv.session.txt`; controls additionally
  `control_indices.txt`.

The engine runs in a unique `.bigwig_attempts/SPACE/METHOD/ID/result` directory.
Only a successful, complete result is renamed into the published method folder.
Each expected file is a Snakemake output, so deleting one track schedules safe
regeneration of that method. Previous publications are archived before execution inside the new
attempt as `previous_publication` (protecting them from failed-job output cleanup); failed attempts retain `driver.log`,
`attempt.json` and engine diagnostics. Do not treat an attempt directory as a
published result. After diagnosing a failure/interruption, resume with:

```bash
capzyme workflow run --route rna --config /data/analysis.yaml --rerun-incomplete
```

The rename publishes within one filesystem. During regeneration the published
method path is absent until success; the former result remains accessible in
the attempt archive. An interruption is recoverable on restart.
Concurrent runs sharing an output directory are unsupported; retain Snakemake's
normal locking. Recovery regenerates the affected method, not individual signal
intervals. Scheduling is one CPU per track job with numerical-library threads
set to one; `mem_mb` is a scheduling request, not an operating-system limit.

Upstream rules no longer depend on every Python script in the repository.
BigWig configuration is resolved separately from upstream job parameters, so
adding methods or changing emitted tracks to compatible completed counts does
not schedule trimming, mapping or counting. Scientific engine/configuration
changes can legitimately invalidate previous jobs; do not disable Snakemake's
integrity/rerun checks to bypass those changes. Older execution snapshots may
have broader recorded dependencies; inspect their preview before upgrading.

See [standalone BigWig semantics](tools/bigwigs.md) for coordinate, control,
float32 and normalization details. Statistical analysis uses the original
integer matrices and has no BigWig dependency.
