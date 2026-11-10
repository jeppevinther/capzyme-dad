# Compact run reports

`capzyme report` summarises an existing analysis directory without rerunning
trimming, mapping, counting or normalization:

```bash
# Regenerate the standard report after inspecting an analysis.
capzyme report --analysis-dir /data/my-analysis/results
# Keep a demonstration or review copy elsewhere.
capzyme report --analysis-dir /data/my-analysis/results --output-dir /data/review/report
```

Default files are `results/report/report.html` and
`results/report/sample_summary.tsv`. `report_status.json` supplies structured
checks and warnings. New workflows also save `analysis_context.json`, the route,
ordered sample IDs, settings and requested-output plan. It records requests,
not execution success. Explicit output directories contain the same report files.
Open `report.html` in a browser. Links are relative to that report directory;
move a report together with its analysis tree or regenerate it after relocation.
A downloaded report cannot open server-only files on a different computer.

The installed command uses Python's standard library only. Reporting requires
neither Snakemake, R, pyBigWig nor MultiQC. Reports may contain private file paths
or sample metadata from saved records; keep them with private analysis results,
not in the public source repository.

## What the report checks

The compact HTML includes run/source/software/reference settings, ordered
per-library accounting, explicit percentage denominators, integrity status,
separate RNA/genome normalization cohorts and factors, and links to matrices,
accepted BAMs, browser tracks and detailed native records. Expand details for
individual checks/output lists. Treatment and replicate fields are included only
when explicitly present in a retained manifest; they are never inferred from IDs.
IDs such as `01` remain strings in the generated TSV (spreadsheet software may
apply its own import formatting).

| Status | Meaning |
| --- | --- |
| Complete and verified (report-time record checks) | A stated small-record checksum, identity, arithmetic, header or presence check passed now; read the stated scope |
| Complete according to recorded checks; no fresh data revalidation | Saved completion/read-identity/input-integrity records report success |
| Not requested / not applicable | Stage does not belong to this route; disabled track methods show “not enabled” |
| Missing, inconsistent or failed | A required receipt/output or consistency check failed |

Overall success means the available saved records pass these checks, not fresh
scientific revalidation. Report generation returns **0** when required checks
pass and **1** when they fail or required outputs are absent. Partial runs still
produce HTML, TSV and JSON explaining failures. Optional unavailable information
is a warning, not failure. File existence alone cannot establish analysis success.

The reporter reads bounded JSON/TSV/text records (16 MiB per record) and matrix
headers. It reconciles exclusive read totals, sample order, QC/publication
bindings, saved handoff pins and overlap results, normalization selection/factor
arithmetic and cross-method source identities. It checks requested output paths.
It does **not** scan FASTQ/BAM bodies, recompute large-file hashes, rejoin read
identities, reconstruct count matrices, decode BigWigs, or invoke scientific tools.
A successful report is therefore not a replacement for an independent validation.
Timing is omitted when reliable start/end records are unavailable; file mtimes
are never used as run times.

## Accounting and denominators

The TSV has one row per sample in recorded order. Blank values mean unavailable
or not applicable, never an inferred zero. Counts are primary input-read outcomes;
nonprimary alignment **records** occupy separate columns and are not added to
input-read totals.

- `raw_reads` and `trimming_removed_reads` come from Cutadapt's input/output counts.
- RNA/genome input, mapped, accepted, low-MAPQ and flag/decoy/CIGAR rejects come
  from stage QC. Final unmapped reads are the last applicable aligner's unmapped
  category; mapped filtering rejects are not silently called unmapped.
- Stage mapped percentage = mapped primary reads / that stage's input reads.
  Stage accepted percentage = accepted reads / that stage's input reads.
  Original-input accepted percentage = accepted reads / raw reads, including
  trimming losses. TSV percentage names spell out each denominator. Undefined
  zero-denominator percentages remain blank.
- In sequential runs, only RNA-aligner-unmapped reads enter genome mapping.
  Raw totals reconcile with trimming losses plus exclusive RNA/genome outcomes.
  Saved overlap=0 and sequence/quality/hash verification are labelled recorded;
  the reporter does not freshly verify read identities.
- Genome-only inputs are already trimmed, so raw/trimming fields are not
  applicable. Existing-BAM inputs have no FASTQ/mapping accounting; saved count
  totals describe the configured alignment strands.
- The older RNA-only streaming route lacks full post-filter category receipts.
  Its native Bowtie2 log supplies mapped/unmapped input accounting; unavailable
  categories remain blank. Both-strand saved count summaries can supply accepted
  totals; one-strand counts must not be labelled complete alignment depths.
- The genome counter's current log has no per-cohort total summary; that limitation
  is stated. Feature values are not rescanned to manufacture a fresh count check.

RPM is position count × 1,000,000 / accepted records in the declared pool.
RNA and genome pools remain separate; one genomic factor applies to both strands.
An explicit external depth pool is labelled as such. DESeq2/control tracks use
count / size factor and state the saved estimator, without invoking R or switching
estimators. Saved method/cohort/units/factor consistency is checked, not reestimated.

## Workflow integration and older runs

All supported routes (`rna`, `sequential`, `genome`, `bam`) now finish with the
same shared report generator. The final rule depends on all requested terminal
analysis outputs, including enabled tracks and optional historical comparison.
No analysis rule depends on reporting. Deleting only a report schedules reporting
without rerunning scientific stages. Adding track methods updates the report while
leaving compatible upstream outputs intact. A no-op workflow restart leaves
reports and analysis outputs unchanged.

The final task cannot certify that Snakemake subsequently exited successfully;
this is explicitly shown. If its required checks fail, Snakemake may remove
failed declared outputs. Diagnostic copies remain in `logs/failed_report.html`,
`logs/failed_sample_summary.tsv` and `logs/failed_report_status.json` beside
`logs/report.log`. Regenerate independently for a partial/incomplete analysis.

Older runs may lack the saved request plan. The reporter infers the route from
saved settings and discovers published normalization records, warning that
unrecorded historical requests cannot be certified. An undiscovered method is
shown as “not enabled in available records (request plan unavailable)”, not a
proven missing requested track. With a request plan, missing enabled methods
fail while disabled methods remain “not enabled.” No manual timing/exit files
are required.

## Optional external MultiQC

MultiQC is not installed, invoked or integrated. Users may run it separately on
recognised native Cutadapt/aligner/SAMtools logs in their own environment.
Keep RNA and genome stages distinguishable (separate stage inputs/reports or
unambiguous stage/sample labels); inspect its recognised modules rather than
assuming every CapZyme accounting record is supported. Keep the original native
logs and this report's explicit accounting alongside any MultiQC summary.

See [implementation and validation results](run_report_validation.md) for the tested
artifact, synthetic routes and read-only completed-run demonstration.
