# Getting started

**Choose a tool:** use the [standalone catalogue](toolbox.md) for individual
operations without Snakemake and the [implemented settings subset](tool_settings.md)
for optional version-2 settings. The tutorial below is the optional coordinator
route, using the unchanged version-1 configuration interface.

Start here with the small **RNA-only synthetic example**. It exercises reference
preparation, trimming, mapping, filtering and 5′ counting without private data,
internet reference downloads, genome indexes, R or institutional access. The
launcher only coordinates existing validation and Snakemake entry points.

## 1. Install and activate the software

Follow the [installation guide](installation.md) to obtain a reviewed checkout
and a working Conda installation. From the source directory:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create --file environment.yaml
conda activate capzyme
python -m pip install --no-build-isolation --no-deps .
capzyme workflow doctor
```

This installs the native tools and Python dependencies together. R is optional
and is not needed for this raw/RPM tutorial. In later sessions, only
`conda activate capzyme` is needed. Platform validation and limitations are in
the installation guide; do not substitute outside tools if resolution fails.

## 2. Create the example and configuration

Run these in the activated terminal. `mktemp` creates a fresh, external working
directory; keep the printed path if you close the terminal. The example refuses
to overwrite a nonempty directory.

```bash
export DEMO="$(mktemp -d "${TMPDIR:-/tmp}/capzyme-tutorial.XXXXXX")"
capzyme workflow example --directory "$DEMO" --with-bigwigs
cd "$DEMO"
pwd -P
capzyme workflow check --config "$DEMO/config.yaml"
```

Expected: `Configuration valid: rna; samples in order: B, A, empty`. The generated
`config.yaml` uses JSON formatting, which is valid YAML. It points to artificial
FASTQs, two tiny FASTAs, `samples.tsv`, and a new `$DEMO/run` output directory.
The third sample deliberately has no retained reads after trimming. Raw/RPM
tracks are enabled for the explicit nonempty B,A cohort; the count matrix still
contains all three libraries. This demonstrates cohort selection without silently
dropping empty samples. Track files are under `run/bigwigs/rna/{raw,rpm}/tracks/`.

For the launcher, the default working directory is the **directory containing
the configuration**. `--workdir /absolute/external/directory` overrides it. Paths
inside the config and sample manifest resolve against that working directory,
not against each individual file's folder. The example writes absolute paths to
remove ambiguity. Direct Snakemake commands instead use their current directory
unless `--directory` is supplied. Keep inputs, outputs and work metadata external
to the source checkout.

## 3. Preview the work

```bash
capzyme workflow preview --config "$DEMO/config.yaml" \
  --cores 4 --memory-mb 16384
```

A preview (Snakemake's **dry run**) lists jobs and dependencies without trimming
or mapping reads. Expect reference-collapse/build, trimming, mapping, counting
and provenance jobs. It may create Snakemake bookkeeping in the external working
directory. Configuration validation checks paths/schema and fixed-reference
checksums where applicable; a preview is not proof that every executable can
process the data. Neither command should create `$DEMO/run/counts/5ends.raw.tsv`.

## 4. Run the example

```bash
capzyme workflow run --config "$DEMO/config.yaml" \
  --cores 4 --memory-mb 16384
```

Wait for successful completion and exit status zero. Snakemake chooses the order
and runs independent jobs together within its scheduling budgets. The values
above are suitable for this tiny example, not a benchmark for a full experiment.
No packages are installed by this command; the activated environment supplies
the tools. See [resources](resources.md) before a substantial run.

## 5. Inspect the result

```bash
head -n 6 "$DEMO/run/counts/5ends.raw.tsv"
cmp "$DEMO/expected.tsv" "$DEMO/run/counts/5ends.raw.tsv"
cat "$DEMO/run/counts/5ends.min_total_10.tsv"
ls "$DEMO/run/qc" "$DEMO/run/logs" "$DEMO/run/provenance"
```

`cmp` prints nothing and exits zero for an exact match. The raw matrix has four
features and three sample columns, in manifest order:

```text
chromosome  strand  start  end  B   A  empty
RNU2__n2    +       1      1    12  0  0
RNU2__n2    +       11     11   0   3  0
virus_minus +      21     21   0   1  0
virus_plus +       10     10   0   1  0
```

The actual file is tab-separated. `start=end` is a one-based nucleotide position.
`virus_minus` is a separate reference identifier; `+` still means forward relative
to that reference. The minimum-total-10 matrix has its header plus only the first
feature (12 counts). The unthresholded matrix remains available. Totals are B=12,
A=5, empty=0. Of 18 raw artificial reads, 17 survive trimming and are accepted;
the adapter-only read is removed. Inspect Cutadapt JSON for input/output counts,
Bowtie2 logs for alignment outcomes and provenance for versions/configuration.

Also expect collapsed reference/membership files in `reference/`, trimmed reads
in `trimmed/`, accepted BAM/BAI files in `bam/`, and `provenance/run.json`. Full
MultiQC aggregation is not implemented. See [results and troubleshooting](results_and_troubleshooting.md)
for route-specific output locations and accounting.

## 6. Resume normally

```bash
capzyme workflow run --config "$DEMO/config.yaml" \
  --cores 4 --memory-mb 16384
```

A completed restart should report **`Nothing to be done`**. Snakemake remembers
which outputs are complete and current. Keep the same external working directory,
its `.snakemake` metadata, configuration and outputs. After an interruption,
diagnose the cause first; then use the same command with `--rerun-incomplete` if
Snakemake reports incomplete files. Do not regenerate the example into its existing
directory to resume. Advanced force-rerun behavior is described separately in
[troubleshooting](results_and_troubleshooting.md#advanced-force-reruns).

## Replace the artificial inputs with your own data

Create a separate external project directory and copy
[core.example.yaml](../config/core.example.yaml) there. Fill all `null` values in
a text editor; do not just change the sample names in the synthetic configuration.
Make a real tab-separated manifest with exactly `sample_id` and `fastq` headers
([template](../examples/core_samples.example.tsv)). Each single-end raw FASTQ is
one library; manifest order becomes matrix-column order.

| Field to review | What you supply/decide |
| --- | --- |
| `samples`, `output_dir` | Manifest path and a new external result directory |
| `reference.ncrna_fasta` | RNA sequences in the intended reference orientations |
| `reference.viral_fastas` | Optional ordered extra FASTAs; use `[]` if absent |
| `prefix_length`, `drop_biotypes` | Default 30-nt collapse and any scientifically justified exclusions |
| `preprocessing` | Your library's adapter, quality threshold, minimum length and NextSeq trim threshold; synthetic settings are examples |
| `mapping.norc` | Whether reverse-complement RNA-reference mappings are forbidden; it does not infer library strandedness |
| `mapping.threads` | Per-mapping CPU request; start with 4 and use a suitable total budget |
| `counting` | MAPQ/duplicate policy, reference-relative strand and optional total-count threshold |

RNA defaults retain MAPQ ≥20, exclude SAM flags 2308, and retain duplicate-marked
and QC-failure reads. Reference collapse groups ncRNAs by their 5′ sequence and
emits a membership table. Review [reference strategy](reference_strategy.md) for
representative selection and identifiers. An existing FASTA/Bowtie2 index instead
uses [fixed.example.yaml](../config/fixed.example.yaml) and checksum pins.

Run `check`, `doctor`, and `preview` against your own configuration before execution.
Obtain the appropriate scientific/resource authorization for experimental runs.
Then choose the relevant route in the [practical usage guide](usage.md).
