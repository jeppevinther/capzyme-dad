# Preserved reference implementations

Three portable originals are included in this repository and remain byte-for-byte unchanged:

- `merge_by_prefix_v2.py`: prefix collapsing, imported by the reference route
  (workflow prefix length 30; original standalone default 70).
- `count5ends_matrix_v1_3.sh`: historical single-nucleotide counting regression
  reference; requires modern Bash and GNU-compatible utilities.
- `make_tss_bigwig_v1_2.sh`: archived visualization reference; its rounding is not
  a required modern normalization mode.

[Public candidate provenance](provenance.public.json) records their SHA-256 hashes.
The project owner confirmed redistribution under the project MIT license on
9 November 2026. The collapse original remains a runtime requirement. See the
[licensing record](../docs/publication/license_review.md).

Other experiment-specific supplied originals and the complete historical
`provenance.json` remain retained for owner review. Private originals are ignored;
do not stage or redact them in place. The public candidate provenance is a
separate subset, not a replacement for the retained record. Nothing in this task
removes previously committed material from Git history.
