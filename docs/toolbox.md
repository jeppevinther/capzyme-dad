# Choose a tool

The scientific tools are the product. Run each independently; Snakemake is an
optional coordinator. Install the [toolbox](installation.md), then set `WORK`
to a separate external project directory. Set `REPO` only for native source-script examples. Activate the dependencies for the selected operation.
No tool installs software or requests experimental execution on your behalf.

| Operation | Program doing the scientific work | Short guide |
| --- | --- | --- |
| Collapse an RNA reference | Preserved `merge_by_prefix_v2.py`, called with a verified source hash | [Reference collapse](tools/reference_collapse.md) |
| Build / verify an index | Native `bowtie2-build` / `hisat2-build` and inspectors | [Index preparation](tools/index_preparation.md) |
| Trim adapters and low-quality ends | Cutadapt | [Trimming](tools/trimming.md) |
| Align to RNA | Bowtie2, with CapZyme filtering, accounting and handoff verification | [RNA alignment](tools/rna_alignment.md) |
| Align to a genome | HISAT2 or Bowtie2, with CapZyme acceptance/QC | [Genome alignment](tools/genome_alignment.md) |
| Count single-nucleotide 5′ ends | SAMtools preparation and BEDTools `genomecov -5 -dz`, shared matrix engine | [FivePrime Count 2](tools/counting.md) |
| Generate strand-specific BigWigs | pyBigWig; optional actual R/DESeq2 factors | [BigWig 2](tools/bigwigs.md) |

```bash
capzyme --help
capzyme count --help
```

Existing script commands remain supported. The optional `toolbox.py` settings
adapter resolves **explicit CLI > version-2 settings > documented defaults**,
prints the exact command, and writes `methods.json` beside each operation's
outputs. Use `--preview` for resolution without execution. The
[implemented settings contract](tool_settings.md) and
[example](../config/tools.v2.example.yaml) describe the supported subset;
the broader version-2 design proposal is not fully implemented.

Mapping adapters do not introduce alignment algorithms. They validate supplied
references/indexes, apply explicit filters/decoy/CIGAR policies, verify exact
aligner-unmapped handoff, and preserve exclusive primary-read accounting. The
native aligners and commands remain visible in `commands.jsonl`. Index creation
is a separate native preparation operation, never an implicit action in this
settings interface. RNA mapping consumes a finished FASTA unchanged; collapse
is a deliberate earlier operation.

Externally generated BAMs are valid counter/BigWig inputs. Dictionaries, actual
sort order, MAPQ, flags, coordinates and source integrity can be checked. Removed
reads cannot be recovered; complete original-read accounting and exact handoff
cannot be reconstructed from an accepted BAM alone. Preserve native logs,
original input identities and all-record alignments if those guarantees matter.

For several tools together, use the existing version-1
[optional workflow launcher](usage.md) and [synthetic tutorial](getting_started.md).
It calls the same implementations. Version-2 files are currently standalone-only;
Snakemake does not silently translate them or introduce a second algorithm.
Optional [workflow BigWigs](workflow_bigwigs.md) reuse unthresholded counts.

Credit upstream software using [software and citations](software_and_citations.md).
Statistical analysis starts from integer counts, never BigWig values.

## Run reporting

`capzyme report --analysis-dir DIR` independently regenerates compact HTML and
TSV summaries from existing records. No scientific stage is rerun. See
[report inputs, outputs, verification scope and examples](run_report.md).
