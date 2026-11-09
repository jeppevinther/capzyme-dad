# Installed toolbox and optional workflow BigWigs

> Historical development snapshot. For the CapZyme-DAD name, current licensing
> and initial-publication checks, see [publication report](publication_report.md).

Snapshot-specific evidence. See [finalization follow-up](finalization_report.md) for
the later dev1 artifact, optional-hash contract and NextSeq-off validation.

6 November 2026. Implementation and synthetic validation only, on
`feature/core-pipeline`, based on `7e13373` plus the reviewed working-tree changes.
No experimental analyses, remote access, research-data changes, commits, pushes,
merges or publication were performed.

## Completed interface

The installable `capzyme-seq` development package supplies `capzyme count`,
`bigwig`, `trim`, `rna-align`, `genome-align`, `collapse`, `index-record` and
`workflow {doctor,example,check,preview,run}`. Existing source-script interfaces
remain supported. Standalone commands do not import or depend on Snakemake.
Installed BigWigs expose matrix, repeated BAM and BAM-manifest inputs; BAM mode
uses the existing validated counter. Native programs/builders remain explicit.

The [installation guide](installation.md) documents base, `tracks`, `settings`
and `workflow` extras. Fresh installation exposed the pyBigWig wheel's NumPy
import requirement; `tracks` now declares NumPy 2 explicitly. Base installation
has no Python runtime dependencies. Simple normalization modes do not invoke R.
A pre-existing adapter bug interpreted `--no-decoys` as false through Python's
negative-boolean option handling. It is now an explicit positive declaration,
covered by a parser regression and actual installed genome mapping with both
aligners. No mapping acceptance or scientific normalization algorithms changed.

Wheel builds include 55 portable runtime files and `source_snapshot.json`, with
per-file SHA-256 hashes. Workflow provenance accepts this installed-source
manifest without inventing a Git identity. Both the wheel manifest and archive
contents were checked against the checkout and for accidental private files.

From a reviewed checkout:

```bash
python3.12 -m venv "$HOME/capzyme-env"
source "$HOME/capzyme-env/bin/activate"
python -m pip install '.[tracks,workflow]'
capzyme workflow example --directory /tmp/capzyme-demo --with-bigwigs
capzyme workflow check --config /tmp/capzyme-demo/config.yaml
capzyme workflow preview --config /tmp/capzyme-demo/config.yaml
capzyme workflow run --config /tmp/capzyme-demo/config.yaml
```

Choose a new/empty example directory and install the documented native programs.
The fixture includes an intentionally empty library: counts preserve B,A,empty,
while its explicit B,A normalization cohort emits four raw/RPM plus tracks.

## Workflow integration

All four routes accept optional version-1 `bigwigs` settings: RNA-only,
sequential RNA-to-genome, genome-only and existing BAMs. The same matrix BigWig
engine receives the raw unthresholded RNA/genome matrix. Thresholded tables are
never substituted. Methods are raw, RPM, DESeq2, controls and supplied factors;
R estimators must be selected explicitly. There is no fallback or inferred control
selection. RNA/genome pools remain separate, with one genomic sample factor for
both strands. Cohort and emitted subset are independent ordered selections.

Automatic RPM denominators audit complete accepted-BAM record depths across all
references and both strands, including when the matrix/display is plus-only.
They do not recount 5′ positions. Explicit audited depth tables are also supported.
Before scheduling, method options, factors, control-table structure, sample
selections, Python dependencies and R/DESeq2 availability are checked. Data-dependent
control matching/informativeness and saved factor/reference/count context are
checked by the original engine once counts are available.

Outputs and exact filenames are documented in [workflow BigWigs](workflow_bigwigs.md)
and [the configuration example](../config/bigwigs.example.yaml). Tracks are under
`bigwigs/SPACE/METHOD/tracks/`, with factors, reference sizes, commands, provenance
and R diagnostics beside them. Unique attempts retain failure evidence under
`.bigwig_attempts/`. Completed method directories are published by rename;
individual tracked files ensure missing outputs schedule safe regeneration.
Old publications are archived before regeneration, protecting them from
Snakemake failed-job output cleanup. The published method path is absent until
success; its former result remains in the attempt archive. An interruption leaves
recoverable attempts; normal locking is required.

Broad upstream code dependencies were narrowed to relevant scientific code;
BigWig configuration is not injected into upstream job parameters. Tests add
tracks to completed compatible RNA and existing-BAM counts and verify upstream
checksums and nanosecond modification times remain unchanged. Older snapshots
may have broader recorded dependencies; inspect their preview on upgrade rather
than suppressing legitimate source-change rerun checks.

## Validation

Complete regression suite: **89 tests in 271.693 seconds; 87 passed, 2 skipped**.
The two skips are the Linux-only process-monitor tests. Fresh-source installed
acceptance tests are enabled within that suite. After the final recovery safeguard, all four workflow/normalization integration
tests passed again in 60.127 seconds, including preservation of a previous
publication on failed replacement. The rebuilt archive-built wheel additionally
passed both installed acceptance tests in 18.868 seconds.
The acceptance coverage includes:

- All four routes; disabled tracks preserve the existing count-only behavior.
- Byte-identical standalone/workflow BigWigs and factors for raw, RPM, actual
  DESeq2-ratio, controls-poscounts and supplied size factors. Existing substantive
  engine tests also cover alternative estimators, sparsity, controls and saved
  factor regeneration.
- Ordered default/explicit cohorts and independent emitted subsets; separate
  RNA/genome scaling; a common factor on both genomic strands; decoded boundary
  coordinates and values; unthresholded versus thresholded source selection.
- Missing dependencies, missing/irrelevant method inputs, invalid factors,
  malformed controls and absent samples fail clearly.
- No-op restart, additional-method generation without upstream changes, a missing
  track, retained failed attempts, a simulated interrupt, and successful recovery.
- Fresh source installation and source-archive-built wheel installation outside
  the checkout. All seven standalone commands execute on tiny synthetic inputs,
  with Snakemake absent from both the standalone environment and PATH. Both
  Bowtie2 and HISAT2 installed mapping adapters pass.
- Exact README two-BAM count and raw-track commands, with temporary paths replacing
  `/data/`, preserve A,B sample order and produce four independently decoded
  expected tracks. All README Bash blocks pass syntax checking. The installed
  synthetic workflow generates the expected count matrix/four tracks and resumes
  unchanged.
- All 99 local links/anchors in the README and updated user-facing guides resolve;
  Markdown fences and navigation ordering pass. GitHub's live rendering has not
  been checked through publication.

### Tested environments

macOS Apple Silicon, Python 3.12.12. Regression coordinator: Snakemake 9.9.0,
PyYAML 6.0.3, pyBigWig 0.3.26. Fresh workflow installation resolved Snakemake
9.27.1, PyYAML 6.0.3, pyBigWig 0.3.26 and NumPy 2.5.4. Both Snakemake versions
passed the relevant workflow tests. Existing native programs: Cutadapt 4.9,
Bowtie2 2.5.4, HISAT2 2.2.1, SAMtools 1.22, BEDTools 2.31.1. Actual R tests used
R 4.5.3 and DESeq2 1.50.2 from the existing isolated local environment.

Fresh standalone environments contained CapZyme, pyBigWig and NumPy, without
Snakemake or PyYAML. Raw/RPM calls explicitly selected a nonexistent Rscript to
verify independence. No native tools, shared environments or R packages were
installed or changed. The installed interface was not executed on Linux in this
task; Linux-only process-monitor tests are skipped on macOS. Native Windows/WSL,
other Python versions and fresh native-tool installation remain untested.

### Review artifacts

Temporary distribution artifacts are under `/tmp/capzyme-final-distributions/`:

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `capzyme_seq-0.1.0.dev0-py3-none-any.whl` | 109827 | `af3969f94ea32803e907e207fc158aec4759a9a68535d785bf53a79b1e9cbe7b` |
| `capzyme_seq-0.1.0.dev0.tar.gz` | 121885 | `17d0a1da7612279d97c60d13f4100b96128e00a10df2debe5e00097f7fd48975` |

Embedded source-manifest SHA-256:
`9b5ec94ec0b7ec3fc88f98a424f121be6cf74e3a1cfa2b306a61970a799fc6e1`.
These are review artifacts, not a release. Temporary installed environments are
`/tmp/capzyme-fresh-base`, `/tmp/capzyme-fresh-full` and `/tmp/capzyme-fresh-wheel`.
Synthetic analysis fixtures are automatically removed; compact local logs and
receipts are retained in ignored `logs/toolbox_interface_2026-11-06/`.

## Preservation and remaining review

`git diff --check` passes. Scientific engines, shared R normalization and archived
legacy scripts are byte-unchanged relative to HEAD. The unrelated Milestone 4A
plan/report match their retained pre-task hashes. The previous README developer
history remains in `docs/development_notes.md`; milestone navigation was removed
from the public front page. Version-1 workflow and optional version-2 standalone
settings remain distinct.

Changed public guides/runtime files pass the sensitive-path/credential scan;
no research FASTQ/BAM/FASTA/index/BigWig is tracked. Existing historical environment
reports retain institutional terminology and are preserved for a separate
publication review. The pending unrelated Milestone 4A report remains untracked;
its historical link is available in the working tree but must be handled explicitly
when preparing a reviewed commit/release. Nothing has been staged or committed.

Conda/Homebrew distribution, a configuration-schema migration, new normalization,
full biological analyses, MultiQC/new downstream statistics and publication are
outside this milestone. Review the working-tree implementation and documentation
before authorizing finalization or release.

## Linux follow-up — 9 November 2026

The reviewed wheel subsequently passed isolated Linux installation and targeted
synthetic workflow/standalone acceptance. See the [Linux acceptance report](toolbox_linux_acceptance_report.md)
for resolved versions, exact numerical/restart checks, resource/storage use and
the known independent-R namespace warning. The 6 November evidence above remains
the historical macOS validation record; scientific implementation is unchanged.
