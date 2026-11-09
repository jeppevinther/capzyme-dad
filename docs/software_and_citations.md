# Software, methods and upstream credit

CapZyme adapters add explicit acceptance rules, single-base matrices, integrity
checks, accounting and reproducible interfaces. Bowtie2 and HISAT2 perform the
alignments; Cutadapt performs trimming; SAMtools handles BAMs; BEDTools computes
aligned 5′-end depth; pyBigWig writes tracks. Cite the programs actually used,
their exact versions, reference/annotation releases and CapZyme source identifier.

| Software | Role | Upstream / citation guidance |
| --- | --- | --- |
| Bowtie2 | Indexing and RNA/unspliced-genome alignment | Langmead & Salzberg (2012), *Fast gapped-read alignment with Bowtie 2*; [official project/citation](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| HISAT2 | Splice-aware genome indexing/alignment | Kim et al. (2019), *Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype*; [official project](https://daehwankimlab.github.io/hisat2/) |
| Cutadapt | Adapter/quality trimming | Martin (2011), *Cutadapt removes adapter sequences from high-throughput sequencing reads*; [official documentation](https://cutadapt.readthedocs.io/en/stable/) |
| SAMtools / HTSlib | Alignment filtering, sorting and validation | Danecek et al. (2021), *Twelve years of SAMtools and BCFtools*; [official publications](https://www.htslib.org/) |
| BEDTools | `genomecov -5 -dz` strand-specific aligned 5′ ends | Quinlan & Hall (2010), *BEDTools: a flexible suite of utilities for comparing genomic features*; [official documentation](https://bedtools.readthedocs.io/en/latest/) |
| pyBigWig / libBigWig | BigWig writing/decoding | Record version and [upstream software](https://github.com/deeptools/pyBigWig), including its license/citation guidance |
| DESeq2 (optional) | Actual R size factors, including `controlGenes` | Love, Huber & Anders (2014), *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*; [Bioconductor package](https://bioconductor.org/packages/release/bioc/html/DESeq2.html); use `citation("DESeq2")` in the installed R environment |
| Snakemake (optional) | Dependency/resource coordination, not a scientific algorithm | [Official citation guidance](https://snakemake.readthedocs.io/en/stable/project_info/citations.html) for the installed version |

For an actual tool run, retain `methods.json` (when using the optional toolbox)
and the engine's `provenance.json`, `commands.jsonl`/`commands.log` or counter
`OUTPUT.jsonl`. These record resolved parameters, exact commands, actual versions,
input/reference hashes, output hashes and status. A command preview is not proof
that an analysis ran. Native-only runs should save command lines, stdout/stderr,
exit statuses, `--version` output and immutable input/index manifests too.

A concise methods record should state: single-end library/trimming settings;
prepared reference identity and annotation provenance; aligner preset/orientation/
seed; MAPQ/excluded flags/duplicate policy; decoy/CIGAR policy; BEDTools aligned
5′ coordinates and conversion; sample order; optional threshold; track method
and denominator/cohort. Record computational settings separately from biology.
External BAMs do not establish the original read accounting by themselves.

BigWig generation from a validated matrix does not remap or recount reads.
Retain the corresponding upstream workflow provenance alongside visualization
records. For current synthetic acceptance evidence see the
[finalization report](finalization_report.md).
