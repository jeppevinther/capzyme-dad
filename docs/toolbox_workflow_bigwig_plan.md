# Installed toolbox and optional workflow tracks

Scope approved 6 November 2026: package existing standalone operations, expose the
workflow launcher, and coordinate the validated BigWig engine from unthresholded
workflow matrices. No experimental execution or scientific algorithm changes.

Implementation:

1. Ship a wheel containing the portable runtime resources and a `capzyme` entry
   point. Keep source-script interfaces. Separate optional track, settings and
   workflow dependencies; native aligners/SAMtools/BEDTools and R remain external.
2. Add a strict, optional workflow `bigwigs` block, independently resolved from
   upstream settings. Validate dependencies and method inputs before scheduling.
   Keep RNA/genome cohorts and emitted samples explicit and separate. Derive RPM
   depths from accepted BAM record totals (never from a strand-limited matrix).
3. Run the existing matrix BigWig engine in unique attempt directories. Preserve
   failed attempts; publish only complete outputs. Track individual files so a
   missing track triggers safe regeneration. Narrow upstream code dependencies.
4. Exercise all normalization modes, every workflow route, installed standalone
   matrix/BAM interfaces, no-op/additional-track restarts, and failure recovery on
   synthetic inputs. Run regression tests and fresh-wheel installation checks.
5. Apply the editorial proposal and align linked user guides. Preserve existing
   user edits and technical history. Report evidence/limitations in `docs/`.

No commit, push, merge or publication is part of this implementation task.

Completed: see [implementation and validation report](toolbox_workflow_bigwig_report.md).
