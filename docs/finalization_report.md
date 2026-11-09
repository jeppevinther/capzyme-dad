# Toolbox finalization — 9 November 2026

> Historical development snapshot. For the CapZyme-DAD name, current licensing
> and initial-publication checks, see [publication report](publication_report.md).

**Passed; stopped before staging or commit.** The approved README and complete
workflow examples are in the authoritative checkout on `feature/core-pipeline`,
based on `7e13373d048ecbd1cd5ae96bf24a3d5d4f209f59` plus reviewed pending work.
No experimental library, shared environment, R package, repository visibility,
remote branch or Git history was changed.

## Implemented behavior and documentation

- Expected reference hashes are optional, binding verification. All observed
  hashes are recorded automatically; they are not trusted historical identity.
  Omitted/null/partial/malformed behavior is specified in
  [reference verification](reference_verification.md). Complete supplied trusted
  manifests remain mandatory to satisfy in full. Fixed RNA still reconstructs
  Bowtie2 sequences and checks IDs/order/content, including same-length mutations.
- `nextseq_trim: null` and `--no-nextseq-trim` omit Cutadapt's extra treatment.
  Ordinary quality trimming remains active. Explicit CLI off/numeric overrides
  settings; omitted standalone settings retain default 20. Numeric zero still
  enables G-tail treatment. Existing numeric behavior is unchanged.
- The approved [README](../README.md) retains the figure, tool order, program
  headings/options tables and standalone-first presentation. The two linked
  [complete examples](../examples/workflows/README.md) enable raw/RPM and show
  optional separate-R DESeq2. RNA-first assignment is optional; genome-only takes
  trimmed input and does not trim. Index building and recording remain separate.
- Conda remains the primary installer. Related tool/config/usage/resource guides,
  tests README, portable project instructions and architecture were made consistent.
  Snapshot-specific reports retain their earlier identities with a follow-up link;
  historical evidence is not rewritten to claim new results. Old instructions are
  retained for review rather than presented as current requirements.

Counting, alignment acceptance, normalization functions, R helper and preserved
original scripts are unchanged. The reference-size helper only changes expected
hash validation; it does not change counting. Snakemake remains an optional
package extra despite inclusion in the recommended convenience environment.

## New package identity

| Item | Identity |
| --- | --- |
| Distribution | `capzyme-seq 0.1.0.dev1` |
| Wheel | `capzyme_seq-0.1.0.dev1-py3-none-any.whl` |
| Wheel SHA-256 | `cfb1aa23706dc96ecd74463b77c5f07308a135a0dceed0c75f724c85175a5b57` |
| `source_snapshot.json` SHA-256 | `539f857d19a4b6c1dd3a968ccce6f3d3102259f54a309c8541faab43e4d8548e` |
| Runtime resource files | 56, each checked against checkout and both installed copies |

This identity replaces the earlier dev0 artifact for these acceptance results;
no old wheel result is attributed to dev1. The working tree is not a new Git
revision. After final edits, every shipped runtime hash still matches this wheel.

Built with the existing complete private Conda environment, then installed into
both existing private platform prefixes without resolving/changing dependencies:

```bash
# From the reviewed checkout; output is outside Git.
python -m pip wheel --no-build-isolation --no-deps --wheel-dir "$DIST_DIR" .
# In each existing private platform environment.
python -m pip install --force-reinstall --no-deps "$WHEEL"
python -m pip check
```

Installed commands ran outside any source checkout. Acceptance used the installed
prefix plus `/usr/bin:/bin`, with user imports disabled and numerical threads one.
All 20 required executable origins belonged to that prefix; native tools and
helpers were owned by the retained Conda package records. Resolved dependency
versions are unchanged from [fresh Conda installation](conda_installation_report.md):
Python 3.12.16, Cutadapt 4.9, Bowtie2 2.5.4, HISAT2 2.2.1, SAMtools 1.22.1,
BEDTools 2.31.1, pyBigWig 0.3.26, NumPy 2.5.4 and Snakemake runtime 9.9.0
(interface-common 1.23.0). The existing Snakemake Conda metadata-version quirk
remains recorded in that earlier report; it was not repaired here.

## Validation

| Check | Result |
| --- | --- |
| Complete automated source suite | 95 tests, 82 passed / 13 explicit skips, 133.075 s |
| Hash/off contract tests | Five passed: omission/null, partial selected pins versus complete manifests, malformed values, mismatches, CLI precedence and real Cutadapt tails |
| Fixed reference suite | Seven passed, including strict historic configs, automatic seven-file provenance, small/large index families, content rejection and immutable inputs |
| Native macOS Apple Silicon installed acceptance | Passed, 71 commands, 53.238 s |
| Linux x86-64 installed acceptance | Passed, 71 commands, 53.206 s |
| Both public examples on each platform | Exact A/B order and integer matrices; sequential 16 raw/RPM strand tracks, direct genome eight; independently decoded positions/dictionaries/float32 values |
| Artificial RNA-only example on each platform | Exact four-feature expected matrix (B,A,empty; totals 12,5,0), four decoded raw/RPM tracks and independent accepted depths |
| Restart and method changes | No-op checksums/mtimes unchanged; added factors preserved upstream and earlier tracks; deleted-track recovery regenerated only affected method and retained previous publication |
| Standalone operations | Counting, matrix/BAM raw/RPM equality, supplied factors without R, tiny Bowtie2/HISAT2 mapping, spliced strands and unmapped handoff passed |
| Integrity | Installed/runtime and transferred inputs match manifests; synthetic originals unchanged; all five supplied originals and unrelated Milestone 4A files unchanged |
| Documentation | Proposed public relative links resolve within candidate tree; no personal review-folder paths; local rendered diagram, headings and tables inspected |
| Git | `git diff --check` passed; index empty; HEAD unchanged |

The public-template fixture has RNA counts A=2/B=1 at position 1 and genome
counts A=B=1 at plus 601/minus 760. RNA RPM is 1,000,000 at its observed end;
genome RPM is 500,000 on each strand using a common both-strand depth of two.
The sequential fixture verifies high-quality terminal G preservation with null
and ordinary removal of four low-quality bases. Adding supplied tracks preserved
75 sequential and 30 direct-genome upstream files. Test harnesses import no
scientific implementation from the checkout.

Explicit suite skips: six optional R/DESeq2 integrations (the supplied-factor
portion passed), three unchanged historical-counter tests needing Bash >=4 on
this Mac, two Linux-only process-monitor tests, and two tests needing separate
base/workflow installation variables. Installed acceptance covered those user
operations separately in the convenience environments, which contain Snakemake;
a new no-Snakemake prefix was not built in this task. No R checks were claimed.
Earlier no-Snakemake and R evidence remains tied to its original artifact.

Linux ran detached at niceness 10 with four CPU affinity cores, a monitored
16 GiB aggregate allowance, 4 GiB synthetic-output ceiling and 10 GiB node-local
free-space floor. Supervision completed exit 0 in 55.278 s with 56 samples,
maximum 1.012 s sampling gap and **315,805,696 bytes (301.176 MiB) sampled peak
aggregate RSS**. Host CPU stayed at/below 1.6%, load at/below 1.64 and available
RAM above 979 GiB. No stop limit was breached. Final fixture directory apparent
size was 40,395,541 bytes; the last periodic sample saw 39,727,827 bytes. These
are sampled/final observations, not continuous peak guarantees. Existing shared
software and research storage/data were not changed or read by the fixtures.

Retained attempts: an old temporary build helper lacked `packaging`; building
with the existing complete environment succeeded without dependency changes.
The first new acceptance fixture had an incorrect synthetic quality-string
length, corrected in the harness before successful Mac/Linux tests. Two added
unit fixtures initially failed on macOS `/var` versus `/private/var` aliases;
resolving their temporary paths fixed those test assumptions. One suite attempt
failed because sandbox access blocked the normal Snakemake cache; the complete
cache-authorized rerun passed. Logs for all attempts remain retained.

## Review inventories and remaining decisions

[Publication review](publication/README.md) provides the exact proposed public
file list, staging list and full per-file disposition, covering all 191 tracked
or nonignored new files. The candidate has **158 files**, including three
conditional legacy originals; **33 retained review records** are withheld.
The proposed local staging list has **88 added/modified files**, including earlier
pending toolbox/Conda work. No files are actually staged or removed. Both unrelated
Milestone 4A documents remain byte-identical and excluded.

The three portable originals need owner [redistribution confirmation](publication/license_review.md).
The collapse original is a runtime dependency, so release is blocked pending that
decision. The current repository/history also contains withheld historical
records: accepting the staging list does not produce a sanitized public history.
Choose and approve a publication/history strategy separately. No deletion or
history rewrite was attempted.

README preview used isolated Markdown-it 4.0.0, Mermaid 10.9.3 and a local
headless browser. The diagram rendered as SVG and eleven tables rendered;
figure, tool table and separate mapping tables were visually inspected. A scrolled
capture was blank; static section captures succeeded. **GitHub-specific rendering
remains pending**, since nothing was pushed/published.

The independent R-checker namespace/stack-imbalance warning remains an explicit
[backlog follow-up](backlog.md); it was not diagnosed or repaired here.

Private artifact/source pins, exact platform paths/commands, versions, failed
attempts, output checksums and Linux supervision receipts remain in ignored
`logs/finalization_2026-11-09/`. `linux_run.private.json` locates the retained remote
inputs/logs and node-local synthetic outputs. Mac outputs and local README
screenshots are external; `access.private.md` records their locations. No generated
fixtures, environment, wheel or sensitive logs belong to the proposed public tree.
