# Project instructions

Follow applicable global instructions as well as these project rules.

## Permanent project organization

- Use this Git checkout as the authoritative location for source code, tests,
  and project documentation. Perform development in this checkout.
- Scientific specifications, implementation plans, and validation reports belong
  in `docs/`. Maintain one authoritative version of each document there.
- Milestone summaries belong in `docs/milestones/` with dated filenames, for
  example `2026-10-29_milestone1.md`.
- Codex session directories and exported reports are temporary snapshots, not
  separate development repositories or authoritative documentation. Do not
  develop in session source copies or create unnecessary duplicate reports.
- Never commit raw sequencing data, BAMs, large outputs, private credentials,
  sensitive logs, or sensitive experiment-specific scripts. Keep analysis data
  and outputs outside the checkout and review all proposed staging.
- Do not perform destructive cleanup of old session files without explicit
  user approval. Preserve those snapshots until cleanup is authorized.

## Scientific development

- Expose important scientific operations as independently executable tools;
  Snakemake coordinates them using the same underlying implementation.

- This repository is intended for public reuse. Keep workflow code independent
  of institutional hosts, usernames, storage paths, credentials, and modules.
- Keep experimental reads, alignments, references, indexes, private sample
  metadata, and large results outside Git and outside the checkout.
- Use experiment configuration for biology and external execution profiles for
  site-specific scheduler settings. Support local execution.
- Preserve supplied legacy scripts unchanged in `legacy_scripts/`; record
  provenance and SHA-256 hashes before developing separate replacements.
- Define and test read orientation, strand, CIGAR handling, coordinate systems,
  multimapping, paired-read counting, and RPM denominators explicitly.
- DESeq2 consumes raw integer counts, never RPM values. Preserve biological
  replicate and treatment metadata and validate the design matrix.
- Prefer small deterministic synthetic tests and pinned per-rule environments.
  Record reference provenance, tool versions, parameters, logs, and exit status.
- Do not submit computational jobs or change existing research data during
  setup or audit. Obtain explicit user authorization before running analyses.
- Do not push, publish, or rewrite legacy scripts without user approval.

## Execution and release boundaries

- Obtain explicit authorization before running experimental analyses, installing
  software, publishing, pushing or merging. Completed validations are not ongoing
  permission to rerun or expand an experiment.
- Use private execution configuration for host selection, resource limits and
  storage locations. Keep site authentication instructions outside public code.
- Follow the selected host's policy. Where shared-server direct execution is
  approved, use the local executor and `nice -n 10`; niceness does not limit RAM.
- Check CPU load, available RAM, disk and applicable quota before launch. Monitor
  the complete process tree and host; stop/report breached limits or sustained
  pressure. Do not change hosts automatically.
- Declare a total CPU and memory scheduling budget with per-rule requests that
  account for aligners, sorting and helpers. Scheduling requests are not OS limits.
- Preserve original inputs/references/indexes. Write to fresh external directories
  and retain commands, versions, checksums, logs, exit status and resource samples.
- Before publication, review the exact proposed file list and Git history for
  private information and confirm legacy-script redistribution permission.
