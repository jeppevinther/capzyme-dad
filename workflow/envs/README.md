Milestone 1 provides pinned top-level per-rule Conda manifests and a coordinator
manifest. Activate the coordinator separately, then optionally use Snakemake's
`--software-deployment-method conda`. No university modules are required.

Local validation used existing CLI tools with the same versions, plus an isolated
Snakemake 9.9.0 environment. Conda deployment and Linux resolution remain untested.
These manifests are not complete transitive/platform locks. Pin containers or
platform-specific environment locks after cross-platform validation.
