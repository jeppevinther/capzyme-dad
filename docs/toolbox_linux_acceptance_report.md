# Linux installed-toolbox acceptance

> Historical development snapshot. For the CapZyme-DAD name, current licensing
> and initial-publication checks, see [publication report](publication_report.md).

Snapshot-specific evidence. See [finalization follow-up](finalization_report.md) for
the later dev1 artifact, optional-hash contract and NextSeq-off validation.

9 November 2026. **Passed numerical installation/workflow acceptance, with a
known warning in the independent R checker.** No experimental input was accessed,
no shared software or scientific implementation changed, and nothing was
committed, pushed, merged or published.

## Reviewed artifact and installation

The original `capzyme_seq-0.1.0.dev0-py3-none-any.whl` (109,827 bytes) was
transferred without rebuilding. Local and remote SHA-256 matched the
[completion report](toolbox_workflow_bigwig_report.md):
`af3969f94ea32803e907e207fc158aec4759a9a68535d785bf53a79b1e9cbe7b`.
Embedded runtime-manifest SHA-256:
`9b5ec94ec0b7ec3fc88f98a424f121be6cf74e3a1cfa2b306a61970a799fc6e1`.
All 55 runtime files matched in both installed environments before and after testing.

Two fresh private virtual environments were created using the existing Python
3.12.16: the wheel with `[tracks,workflow]`, and separately `[tracks]` only.
Both passed `pip check`; resolved package reports and complete freezes were
retained. No source checkout was used for execution.

| Component | Resolved version |
| --- | --- |
| Platform | Red Hat Enterprise Linux 8.10, x86-64 |
| Python / pip | 3.12.16 / 25.0.1 |
| CapZyme / Snakemake | 0.1.0.dev0 / 9.27.1 |
| pyBigWig / NumPy / PyYAML | 0.3.26 / 2.5.4 / 6.0.3 |
| Existing Cutadapt / Bowtie2 | 4.4 / 2.4.2 |
| Existing SAMtools / BEDTools | 1.21 / 2.31.0 |
| Existing isolated R / DESeq2 | 4.5.3 / 1.50.2 |

Actual executable paths, native modules, package locations and checksums are in
private receipts. Only the new environments were installed; no R packages or
reused software prefixes were installed into or changed.

## Acceptance results

All **48 recorded test commands exited zero**; installation and the monitored
launcher also exited zero. The portable [checker](../tests/linux_installed_acceptance.py)
uses supplied Linux-native paths without macOS installation assumptions.

| Check | Result |
| --- | --- |
| Installed top-level and all seven individual tool help interfaces | Passed in both environments |
| `capzyme workflow doctor` | All selected dependencies passed |
| `example --with-bigwigs`, configuration check, preview, execution | Passed; preview created no outputs |
| Exact count matrix versus fixture and literal independent expectation | 4 features × 3 samples, all integer cells exact |
| Four expected raw/RPM tracks | Exact dictionaries/order, coordinates and 8 nonzero float32 values; no extra intervals |
| Second normal execution | No jobs; all 67 result-file checksums/sizes/nanosecond mtimes unchanged |
| Add supplied-factor tracks to completed counts | 47 upstream files and existing raw/RPM files unchanged |
| Delete one synthetic raw track and restart | Only affected method regenerated; exact track restored; prior publication archived; other outputs unchanged |
| Standalone count and matrix/BAM BigWigs without Snakemake | Four matrix/BAM track pairs byte-identical and independently decoded; source BAM/FASTA unchanged, no input indexes created |
| DESeq2-ratio and controls-poscounts workflow integration | Native R factor agreement within 1e-12; standalone/workflow tracks byte-identical and independently decoded |

The main matrix preserves sample order **B, A, empty**, with totals **12, 5, 0**.
Its SHA-256 is
`792cd8a5d0bf3debe385035c42f2ece7a2285a52a2f012497a3ce0bf56dae001`.
The example explicitly selects B,A as its normalization/emission cohort; the empty
library remains in the count matrix. Independent accepted-BAM audits returned
B=12 and A=5 across the complete reference/strand pool.

| Sample | Reference | 1-based position, strand | Raw | RPM |
| --- | --- | --- | ---: | ---: |
| B | RNU2__n2 | 1, + | 12 | 1,000,000 |
| A | RNU2__n2 | 11, + | 3 | 600,000 |
| A | virus_minus | 21, + | 1 | 200,000 |
| A | virus_plus | 10, + | 1 | 200,000 |

Each BigWig interval is `[position-1, position)`. Strand is relative to the
reference record, including the synthetic `virus_minus` record. The separate
standalone fixture verifies both strands at positions 1 and 100, with raw value 1
and RPM 500,000 from two accepted alignments.

R integration uses a separate two-library/two-position synthetic BAM fixture,
with cohort B,A but tracks emitted only for A on both strands. Both methods
agree with independent **native DESeq2** calls. Controls use two explicitly
selected synthetic positions and minimum 2 for this fixture; production defaults
are unchanged. Thresholded outputs are not used for normalization.

## R warning and coverage limits

Both independent native R checks reproduced the previously reviewed first-use
namespace warnings: stack imbalance in `::` and `<-`. All numerical comparisons
passed. The four production normalization command logs—workflow and standalone,
both methods—contained no R warnings. This matches the earlier cold-namespace
observation; it is not evidence that the R environment is universally warning-free.
The underlying initialization issue remains unresolved. No additional diagnostic
analysis, package repair, estimator fallback or version substitution was attempted.

There was **no unavailable-R coverage skip**. This targeted Linux acceptance
covers the RNA example and existing-BAM normalization route. Full Linux regression,
genome-only/sequential BigWig orchestration, other estimator/control edge cases,
other Python/platform combinations and experimental execution were not repeated.

## Resources, integrity and retained evidence

Detached execution ran **10:15:13–10:24:30 CET** on the approved shared compute
host, with `nice -n 10`, four-CPU affinity, local executor, `mem_mb=16384`, a 16 GiB
process-tree RSS guard and one numerical-library thread. Niceness is not a memory
limit. Installation plus acceptance took **555.374 s** under monitoring; the
acceptance driver itself took **245.910 s**.

There were **556 resource samples**, maximum wall-clock gap **1.022 s**, and
sampled peak aggregate RSS **1.542 GiB**. No resource stop condition occurred.
Host CPU reached 2.1%, load reached 10.97, and available RAM stayed above
982 GiB. These are sampled observations, not continuous peak guarantees.

Retained project-file apparent size before the final receipt was **283,755,380
bytes (270.610 MiB)**, including both environments and all synthetic outputs/logs.
The storage guard ceiling was 20 GiB with a 50 GiB filesystem-free floor; no space
reservation was allocated. Project quota capacity was unobservable (`quota -s`
exited 1 without capacity information); filesystem free space is not quota
confirmation. Node-local pip scratch storage was not continuously measured.

Wheel/checker integrity, both installed runtime manifests, synthetic input
integrity, and eight reused executable checksums passed final verification. The
detached session has ended. Initial authentication and read-only preflight logging
needed retries; these preceded installation. The Linux installation/acceptance
run itself required no recovery beyond the deliberately tested missing track.

Complete commands, stdout/stderr, statuses, software paths/freezes, source pins,
monitor samples and synthetic outputs remain in the new private Linux directory.
The exact location is in ignored
`logs/linux_toolbox_acceptance_2026-11-09/access.private.md`; compact local receipts
are in its `linux_receipts/` directory. Only compact logs and summaries were
transferred back (273,327-byte archive plus four small normalization logs), not
FASTQs, BAMs, matrices, BigWigs or environments. Git staging remains empty.
