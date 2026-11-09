# Synthetic regression tests

From the reviewed checkout, with the recommended environment activated:

```bash
python -m unittest discover -s tests -v
```

Fixtures are generated in external temporary directories; no experimental inputs
are needed. Native integrations require Cutadapt, Bowtie2/HISAT2 and their builders/
inspectors, SAMtools, BEDTools, pyBigWig and Snakemake. Missing optional tools cause
explicit skips. Set `CAPZYME_TEST_BASH` to Bash >=4 for the unchanged historical
shell-counter regression (GNU-compatible sort/awk also required). Linux process-tree
monitor tests require Linux/psutil. Optional actual R tests use
`CAPZYME_BIGWIG_RSCRIPT` with a compatible DESeq2 environment.

`test_fixed_reference.py` covers strict/optional hashes, read-only reuse and strong
sequence compatibility. `test_finalization.py` covers hash contracts and the true
NextSeq-off mode, including CLI precedence and real Cutadapt tail trimming.
Existing tests cover strands/coordinates/filtering, mapping/splicing/decoys,
handoff provenance, accounting, track normalization and safe restart/recovery.

The installed acceptance drivers import no scientific code from the source tree.
`finalization_acceptance.py` exercises the actual wheel outside any checkout,
including both complete public workflow templates, raw/RPM tracks, restart and
method additions. It requires a reviewed source-manifest SHA-256 argument.
`conda_install_acceptance.py` additionally verifies all native executable origins.
See the [finalization report](../docs/finalization_report.md) for exact artifact
identities, commands, passed cases and untested optional coverage. Historical
milestone test totals describe their earlier snapshots, not the current suite.
