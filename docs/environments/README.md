# Tested Conda environment records

Use the readable [root environment recipe](../../environment.yaml) for the
recommended installation, followed by the explicit source-install command in the
[installation guide](../installation.md). These platform-specific explicit Conda
records capture exact tested builds/URLs; they complement that recipe rather than
replacing it with an opaque machine export. No user prefixes or credentials are
included.

- `conda-osx-arm64.explicit.txt`: native Apple Silicon package record.
- `conda-linux-64.explicit.txt`: Linux x86-64 package record.

CapZyme itself is installed separately from reviewed source, so it is not a
Conda channel entry in these files. Exact source/runtime identities, executable
origins, tests and limitations are in the [acceptance report](../conda_installation_report.md).
Do not use another architecture's record or substitute emulation.
