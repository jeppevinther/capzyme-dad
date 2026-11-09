# Installation

The recommended environment installs the native programs and Python dependencies
together. R/DESeq2 remains optional. This installs CapZyme-DAD from source into your
own environment; it is not a published Conda/Bioconda package.

## 1. Have Conda available

Use an existing, initialized Conda installation. If you need one, follow the
[official Miniforge installer instructions](https://github.com/conda-forge/miniforge#install),
choosing your actual operating system and CPU architecture. Open its initialized
terminal before continuing. The installer and its base environment are separate
from the CapZyme environment below.

## 2. Download CapZyme-DAD

Download and unzip the [source archive](https://github.com/jeppevinther/capzyme-dad/archive/refs/heads/main.zip), or use Git:

```bash
# Download the software source into a new folder.
git clone https://github.com/jeppevinther/capzyme-dad.git

# Enter the folder containing environment.yaml and pyproject.toml.
cd capzyme-dad
```

If you downloaded the ZIP instead, open a terminal in the extracted folder.
CapZyme-DAD installs from these source files. Its installed command is `capzyme`.

## 3. Create and activate the environment

From the source directory, run:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create --file environment.yaml
conda activate capzyme
```

The environment file uses conda-forge, then Bioconda, with default channels
excluded. Strict priority applies to this creation command without changing your
global Conda configuration. Choose a new environment; do not install into `base`
or overwrite an existing `capzyme` environment.

## 4. Install CapZyme-DAD

With `capzyme` activated and still in the source directory:

```bash
python -m pip install --no-build-isolation --no-deps .
```

Conda has already installed the build requirements and dependencies. These flags
use those build tools and prevent pip from replacing Conda's packages. This is a
regular installation, not an editable/source-path installation. Installed commands
work from any directory. The final explicit step avoids ambiguity about where a
relative pip entry inside an environment file would run.

In later terminal sessions, simply run `conda activate capzyme`. You do not need
to recreate the environment or reinstall CapZyme. Reinstallation is needed only
when deliberately installing a different reviewed source version.

## 5. Check installation

```bash
capzyme --help
capzyme workflow doctor
capzyme workflow doctor --route genome --aligner bowtie2
capzyme workflow doctor --route genome --aligner hisat2
python -m pip check
```

These commands check the installed interfaces and dependencies; they do not run
an experiment or create a genome index. Builders and inspectors for both aligners
are included. The launcher uses the activated environment; no per-rule environment creation
flag is needed. The convenience environment contains Snakemake, but standalone
CapZyme commands do not require or invoke it.

## 6. Run the small synthetic example

Choose a new or empty directory outside the source checkout:

```bash
capzyme workflow example --directory /tmp/capzyme-demo --with-bigwigs
capzyme workflow check --config /tmp/capzyme-demo/config.yaml
capzyme workflow preview --config /tmp/capzyme-demo/config.yaml --cores 4 --memory-mb 16384
capzyme workflow run --config /tmp/capzyme-demo/config.yaml --cores 4 --memory-mb 16384
```

The artificial example produces:

- Counts: `/tmp/capzyme-demo/run/counts/5ends.raw.tsv` (B, A, empty; totals 12, 5, 0).
- Expected counts: `/tmp/capzyme-demo/expected.tsv`.
- Four tracks: `run/bigwigs/rna/raw/tracks/{B,A}.plus.bw` and
  `run/bigwigs/rna/rpm/tracks/{B,A}.plus.bw` under that example directory.
- Factors, reference sizes, commands and provenance beside each method's tracks.

The example explicitly uses the two nonempty libraries as its normalization
cohort. Repeating the `run` command should report nothing to do. No R, research
libraries or reference downloads are needed. See the [tutorial](getting_started.md)
for count coordinates and [workflow tracks](workflow_bigwigs.md) for other methods.

## Platforms and reproducibility

Tested platforms: **Linux x86-64 (RHEL 8.10)** and **native macOS Apple Silicon**.
Windows/WSL, macOS Intel and Linux ARM remain untested.

Platform acceptance and resolved package records are documented in the
[Conda installation report](conda_installation_report.md) and
[platform environment records](environments/README.md). Support statements follow
actual installation/execution results; a solver success alone is insufficient.
The readable [environment.yaml](../environment.yaml) is the primary recipe;
platform-specific records complement it. Do not use architecture emulation or
outside tools to conceal an unavailable dependency.

## Optional R/DESeq2

| BigWig normalization | R needed? |
| --- | --- |
| `raw`, `rpm`, `supplied` | No |
| `deseq2`, `controls` | Yes: a separately installed compatible R/Bioconductor DESeq2 environment |

Use an explicit `--rscript` for the optional R environment. It is not installed,
repaired or upgraded by the recommended recipe. Statistics use the original
integer count matrix, never BigWig values. See [BigWig normalization](tools/bigwigs.md).

## Advanced: minimal pip installation

For users who already manage the native software, Python ≥3.12 can install
CapZyme in a separate virtual environment. This alternative does not install
native aligners, SAMtools or BEDTools:

```bash
python3.12 -m venv /path/to/new-capzyme-venv
source /path/to/new-capzyme-venv/bin/activate
python -m pip install '.[tracks,workflow]'
```

Run the final command from the reviewed checkout. Optional groups are `tracks`
(pyBigWig 0.3.26 and NumPy 2), `settings` (PyYAML), and `workflow` (Snakemake and
PyYAML). `pip install .` has no Python runtime dependencies; standalone use does
not acquire Snakemake as a mandatory dependency. Existing script interfaces remain
supported. Do not mix this dependency-resolving alternative into the recommended
Conda installation in place of its explicit `--no-deps` step.

| Operation | Software for independently managed installations |
| --- | --- |
| Collapse | Python; preserved collapse code ships with CapZyme |
| Index building/recording | Native Bowtie2/HISAT2 builder and inspector |
| Trimming | Cutadapt |
| RNA/genome alignment | Bowtie2 or HISAT2, plus SAMtools |
| Counting / BAM BigWigs | SAMtools and BEDTools |
| Matrix BigWigs | pyBigWig and NumPy; R/DESeq2 only for the two R modes |

Please credit [the upstream software](software_and_citations.md) used in your analysis.
