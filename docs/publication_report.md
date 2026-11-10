# CapZyme-DAD initial public source

The project owner selected **CapZyme-DAD — Data Analysis and Display for
CapZyme-seq** and authorized publication on 9 November 2026. The installable
package is `capzyme-dad`; the installed command and Conda environment remain
`capzyme`, preserving the documented user interface.

## Public tree and licensing

The public repository is [jeppevinther/capzyme-dad](https://github.com/jeppevinther/capzyme-dad).
It starts from a reviewed portable source export with a fresh Git history.
The internal development repository, research data, credentials, private logs
and withheld historical records were not copied into that history.
The initial export contained 159 files. The
[public source inventory](publication/public_files.txt) tracks the current tree,
including additions made after that initial publication.
The owner confirmed redistribution of the three original portable legacy scripts
under the project MIT license; their hashes and author notices are unchanged.
See [licensing](publication/license_review.md).

## Renamed package identity

| Item | Value |
| --- | --- |
| Distribution | `capzyme-dad 0.1.0.dev2` |
| Wheel | `capzyme_dad-0.1.0.dev2-py3-none-any.whl` |
| Wheel SHA-256 | `03d6659c2ce81002fa7463032e658cb0a3bdffab6e91714773710e43681466a0` |
| Runtime manifest SHA-256 | `43350d29f5d21f261e6704ad28ed420556ac0f72c9461ecfb783e2ce613cf6a2` |
| Hashed runtime resources | 56 |

The package and its embedded source manifest use the new name. Earlier dev0/dev1
acceptance reports retain their original artifact identities and are not relabeled
as tests of this wheel. This is a public source publication, not a PyPI/Bioconda
submission or a claim of a stable production release.

## Verification before upload

- The source suite ran **95 tests: 82 passed, 13 explicit skips**, in 145.169 s.
  Skips match unavailable optional R/historical Bash/platform-specific or separate
  installation fixtures; no skipped coverage is reported as passing.
- An isolated regular installation of the renamed wheel passed **71 synthetic
  acceptance commands**, all exit zero, in 55.815 s. It ran outside
  the checkout and reused the previously validated private macOS Apple Silicon
  Conda dependency stack. No fresh dependency installation is claimed here.
- Both complete public examples passed: sequential RNA-to-genome mapping and
  direct-genome mapping, exact count matrices, decoded raw/RPM strand tracks,
  NextSeq-off behavior, restart preservation, extra methods and recovery.
- All 56 installed runtime hashes matched the embedded manifest. Wheel and
  source identities were recorded separately from prior validation artifacts.
- The preserved original scripts match their recorded SHA-256 hashes. Local
  documentation links resolve; no personal paths or credential patterns were
  found in the export. Original/private Git history was not copied.

Counting, alignment acceptance and statistical normalization algorithms were not
changed by the naming/publication step. The scientific implementations and
Conda dependency versions retain their earlier Linux/macOS evidence, including
the [finalization tests](finalization_report.md). This renamed wheel was not
independently reinstalled on Linux during publication; its metadata/helper
acceptance checks were updated for the new distribution name.

The independent R-checker warning remains a [documented follow-up](backlog.md).
No R environment or experimental library was changed or rerun. Build artifacts,
private command receipts and synthetic outputs are kept outside this repository.
