# Reference identity: observed hashes and optional expectations

Ordinary use requires existing reference/index paths and substantive compatibility
checks. It does not require copying hashes manually. Observed SHA-256 values are
recorded automatically; observing a hash now does not prove historical identity.

For fixed RNA workflow references, `fasta_sha256` and `index_sha256` are optional.
Omitted or null fields mean no expected pin; `{}` also means no index expectations.
A partial index map is allowed: every supplied key must be one of the selected
six `.bt2`/`.bt2l` suffixes and every supplied value must be 64 hexadecimal digits.
An individual null map entry is invalid; omit an unpinned key instead. Unknown
suffixes, non-map index expectations, empty/invalid hashes, and non-string values
fail before jobs are scheduled. Upper/lowercase hex is accepted. Every supplied
expectation is binding; a mismatch never becomes a newly trusted observation.

Fixed RNA verification still reconstructs Bowtie2 sequences and compares reference
identifiers, order and sequence content (IUPAC ambiguity masked to N only for this
comparison). It hashes all seven input files before/after and records the full
observed map plus supplied expectations in `reference/fixed_verified.json`.
Supplied references/indexes are never regenerated or modified.

Optional strict example:

```yaml
reference:
  mode: fixed
  fasta: /data/rna.fa
  index_prefix: /data/rna-index
  index_format: bt2
  fasta_sha256: REPLACE_WITH_64_HEX_DIGITS
  index_sha256:
    .1.bt2: REPLACE_WITH_64_HEX_DIGITS
```

This example pins only the FASTA and first index part. Supply the remaining
suffixes to pin the complete family. Obtain trusted pins from retained earlier
build records, not by accepting a changed current file.

Existing-BAM `reference.fasta_sha256` and genome `reference_sha256` follow the same
omitted/null/valid-string contract. Their observed reference hashes are retained
in run/stage provenance. Historical matrix comparison retains its independent
strict expected-matrix checksum contract.

Standalone mapping settings `references.SPACE.index.verification_manifest` and
workflow genome `verification_manifest` are optional (omitted/null = no trusted
manifest). Without one, mapping still validates the FASTA/index dictionary and
records observed FASTA/index hashes and read-only checks in native provenance.
`capzyme index-record` can record an existing index without claiming trusted
historical origin. With `--trusted-manifest` it must match the earlier record.

A supplied trusted manifest must be complete: schema 1, correct aligner, a valid
FASTA hash and all six/eight index-part hashes. Partial, null or malformed trusted
records fail; this differs intentionally from selected individual expected pins.
Name/length compatibility is not proof of sequence or embedded annotation origin;
retain native build/annotation records. Fixed RNA's stronger content check is
preserved separately. See [index preparation](tools/index_preparation.md).
