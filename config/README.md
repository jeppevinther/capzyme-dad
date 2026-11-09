# Tools, workflows and examples

For complete editable workflow settings and sample sheets, start with
[the two public examples](../examples/workflows/README.md): sequential RNA→HISAT2
from raw reads, or direct Bowtie2 genome mapping from trimmed reads. Both enable
raw/RPM tracks and show optional DESeq2 tracks with a separate R environment.

| File | Current purpose |
| --- | --- |
| `core.example.yaml` | Version-1 generated RNA reference workflow |
| `fixed.example.yaml` | Version-1 existing RNA FASTA/Bowtie2 index; no required expected hashes |
| `genome.example.yaml` | Version-1 genome-only route; already trimmed FASTQs |
| `existing_bam.example.yaml` | Version-1 existing-BAM counting; no trimming/mapping |
| `bigwigs.example.yaml` | Optional version-1 workflow track block |
| `tools.v2.example.yaml`, `tool_samples.example.tsv` | Optional standalone version-2 settings/samples; not workflow configs |
| `bigwig.deseq2.v2.example.yaml` | Optional standalone DESeq2 track settings |
| `config.example.yaml` | Historical broad design sketch; not accepted by any current workflow |

Copy templates to an external analysis folder and replace required input/output
paths. In the workflow schema, relative paths resolve against the working directory
(the configuration directory with the launcher). In standalone settings, paths
resolve against the settings file; sample input paths resolve against its sheet.
The two schemas are intentionally distinct.

Observed reference/index hashes are automatic provenance. Expected FASTA hashes
and selected fixed-index pins are optional, binding verification; missing/null
fields mean unpinned. Null values for required file locations remain invalid.
See [the exact hash contract](../docs/reference_verification.md), including partial
expectations and complete trusted manifests. Content/dictionary checks remain
active without expected hashes; supplied indexes are never rebuilt implicitly.

`preprocessing.nextseq_trim: null` omits the extra two-color Cutadapt treatment
while ordinary quality trimming remains active. Numeric zero still enables it.
Choose a numeric threshold for suitable two-color NextSeq/NovaSeq chemistry.

Historical matrix comparison is a separate optional exact-comparison contract,
including a trusted expected-matrix checksum; do not copy an experiment's hashes
into generic reference settings. Keep real inputs, private paths and outputs
outside Git and outside the source checkout.
