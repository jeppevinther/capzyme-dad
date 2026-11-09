# Analysis examples

Use the [complete workflow settings and sample sheets](workflows/README.md) for
sequential RNA→HISAT2 mapping from raw reads or direct Bowtie2 genome mapping
from trimmed reads. Both enable raw/RPM tracks. Replace their placeholders in
an external analysis folder before running; no experimental data are supplied.
For a ready-to-run artificial RNA example use `capzyme workflow example
--directory /tmp/capzyme-demo --with-bigwigs` in a new/empty external directory.

The older `samples.example.tsv` is a future-design sketch, not accepted workflow
input or a valid differential design. `core_samples.example.tsv` illustrates the
implemented two-column single-end sheet. `rsv_parameters.yaml` is an explicit
protocol preset, not real sample metadata or a universal trimming recommendation.
Keep real experiments and biological replicate/treatment metadata outside Git.
