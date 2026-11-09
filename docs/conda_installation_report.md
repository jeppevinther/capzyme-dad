# Complete Conda installation acceptance

> Historical development snapshot. For the CapZyme-DAD name, current licensing
> and initial-publication checks, see [publication report](publication_report.md).

Snapshot-specific evidence. See [finalization follow-up](finalization_report.md) for
the later dev1 artifact, optional-hash contract and NextSeq-off validation.

9 November 2026. **Passed on Linux x86-64 and native macOS Apple Silicon.**
The recommended route now installs native tools, Python dependencies and build
requirements together, followed by a regular CapZyme source installation.
No scientific implementation, interface or R environment changed. No experimental
data, package submission, publication, commit or shared-environment installation
was involved.

## Exact recommended commands

From the maintainer's reviewed source directory, with Conda already initialized:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create --file environment.yaml
conda activate capzyme
python -m pip install --no-build-isolation --no-deps .
```

The explicit pip step was chosen over a YAML pip section to make its source
working directory unambiguous and prevent dependency replacement. Conda supplies
Python, pip, setuptools and wheel. The package remains non-editable and works
outside the checkout. Snakemake is included for convenience but remains an
optional package dependency; `pyproject.toml` is unchanged.

Both tests used clean new prefixes and private package caches; CONDA_ENVS_PATH
redirected the environment name into private test storage without changing the
documented creation/activation/install commands. The Mac final
prefix reused only downloads from this task's initial attempt, not an installed
environment or its tools. Linux installed the same source/build files from a
hash-verified 61-file temporary deployment snapshot. No development occurred in
that snapshot. Scientific runtime manifest SHA-256 remained
`9b5ec94ec0b7ec3fc88f98a424f121be6cf74e3a1cfa2b306a61970a799fc6e1`.

## Resolved software and origins

| Component | Linux x86-64 | macOS Apple Silicon |
| --- | --- | --- |
| Conda used to create environment | 24.5.0 | 26.9.4 |
| Python | 3.12.16 | 3.12.16 |
| Cutadapt / Bowtie2 / HISAT2 | 4.9 / 2.5.4 / 2.2.1 | 4.9 / 2.5.4 / 2.2.1 |
| SAMtools / BEDTools | 1.22.1 / 2.31.1 | 1.22.1 / 2.31.1 |
| pyBigWig / NumPy | 0.3.26 / 2.5.4 | 0.3.26 / 2.5.4 |
| Snakemake runtime / interface-common | 9.9.0 / 1.23.0 | 9.9.0 / 1.23.0 |
| setuptools / wheel | 84.2.1 / 0.48.0 | 84.2.1 / 0.48.0 |

Strict priority was used with conda-forge first, Bioconda second and default
channels excluded. All installed Conda records came from those two channels.
Linux Conda printed inherited default-channel notices, but installed no default
channel packages. Records contain only `linux-64`/`noarch` or
`osx-arm64`/`noarch`, respectively; no emulation was used.

All 20 checked executables resolved within the new prefix. Native programs,
small-index builders/inspectors, annotation extractors, Python and Snakemake
were owned by that environment's Conda records; CapZyme was pip-installed there.
Acceptance PATH was restricted to the new prefix plus `/usr/bin:/bin`. Existing
server modules, system bioinformatics tools and Homebrew could not mask missing
programs. HISAT2 inspection uses its native inspector binary through the existing
validated adapter; no vendor wrapper was rewritten.

Exact build/URL records are in [docs/environments](environments/README.md), separately
for each platform. Private receipts retain full executable paths, resolved
package JSON, source pins and installation logs.

## Synthetic acceptance

Both platforms passed **54 recorded commands**, all exiting zero:

- Top-level/individual installed help and dependency doctors for RNA, Bowtie2
  genome and HISAT2 genome routes.
- The installed `example --with-bigwigs`, configuration check, preview and run.
  The 4-feature × 3-sample matrix was byte-identical to the fixture: ordered
  samples B,A,empty; totals 12,5,0; SHA-256
  `792cd8a5d0bf3debe385035c42f2ece7a2285a52a2f012497a3ce0bf56dae001`.
- Four raw/RPM tracks with exact independent reference dictionaries, single-base
  coordinates and all eight nonzero float32 values. A second execution scheduled
  no work and preserved all 67 output hashes/sizes/nanosecond mtimes.
- Direct standalone counting, raw/RPM matrix-versus-BAM track agreement, both
  strands and supplied-factor tracks, with a deliberately unavailable Rscript.
- Tiny native Bowtie2 mapping: two accepted reads, one unmapped; exact 5′ positions
  `chr:+:601` and `chr:-:760`.
- Tiny native HISAT2 mapping with generated annotation/index: four accepted reads,
  one unmapped; the above ends plus `chr:+:171`, `chr:-:470`; both splice CIGARs
  exactly `50M200N50M`. Native builders, index recording/inspection, counting and
  unmapped handoff were exercised. Synthetic source inputs/indexes were unchanged.
- Additional-track isolation and missing-track recovery also passed through the
  reused installed-interface checker.

The [Conda checker](../tests/conda_install_acceptance.py) reuses the earlier
[installed acceptance checker](../tests/linux_installed_acceptance.py). It imports
no source-tree scientific implementation. R/DESeq2 was deliberately absent and
not tested or repaired; raw/RPM/supplied modes passed without it.

Focused regressions in the final Mac Conda environment: **40 tests in 57.478 s;
36 passed, 4 R-dependent tests skipped explicitly**. These cover toolbox,
launcher, genome mapping and BigWig behavior. Complete old scientific validation
and experimental-library execution were not repeated.

## Failures, resources and limitations

The initial Mac recipe solved and passed synthetic execution but **failed
`pip check`**: interface-common 1.23.1's Python metadata requires packaging ≥26.1,
while Snakemake 9.9's Conda constraints require packaging <26. The Conda recipe
for interface-common did not express that higher requirement. The final readable
recipe explicitly pins the compatible 1.23.0; a second clean prefix passed
installation, `pip check` and acceptance. The failed first environment/logs remain
untouched. No third-party metadata or scientific behavior was patched.

The Bioconda Snakemake package reports runtime/Conda version 9.9.0 but Python
package metadata version 0.0.0. Runtime/version checks and `pip check` passed.
The documented `--no-deps` install avoids pip replacing Conda's selected stack;
advanced pip-extra installation is a separate route, not an upgrade instruction
for this environment.

An early Linux transfer raced directory creation. The checksum gate prevented
any installation; transfer was retried after creation and verified. The final
Linux installation/acceptance run exited zero without runtime recovery.

Linux ran detached under nice 10, four-CPU affinity, a 16 GiB process-tree RSS
guard and one numerical-library thread. Full solve/install/acceptance took
116.179 s; acceptance itself 21.918 s. Sampled peak aggregate RSS was 1.033 GiB;
117 samples had maximum gap 1.025 s and no stop conditions. Native Apple Silicon
acceptance took 33.025 s. Sampling does not prove continuous peak usage.

Linux used new private node-local prefixes/cache for bounded metadata I/O, with
durable logs in private project storage. Sampled logical file size including
cache/prefix was 2.076 GiB (hardlinked files can be counted twice), below a 10 GiB
ceiling; the filesystem-free floor was 10 GiB. It was not a disk reservation or
project-quota confirmation. Existing research data/software were untouched.

Windows/WSL, macOS Intel, Linux ARM and other package-version combinations are
untested. The readable recipe constrains direct versions; exact platform records
retain the tested transitive builds. Re-solving in future can change transitive
packages and permitted patch versions. Nothing was published to Conda/Bioconda or a package index.

## Files and retained records

Changed: root `environment.yaml`; README installation section;
`docs/installation.md`; `docs/getting_started.md`; this report;
`docs/environments/` platform records; and the new Conda acceptance checker.
Scientific source, package dependency declarations and existing interfaces are
unchanged. Unrelated pending work is preserved and staging is empty.

Commands, statuses, origin hashes, full package lists and logs are retained in
ignored `logs/conda_acceptance_2026-11-09/`; exact environment/output locations are
in its private access receipt. Environments and generated synthetic data remain
outside Git. Stop for review before any finalization/publication.
