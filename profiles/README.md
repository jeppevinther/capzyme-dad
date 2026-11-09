The local profile supports scaffold dry runs. The SLURM profile requires a
compatible executor plugin, scheduler commands, an authorized submission host,
and an external user profile with site-specific settings. Its resource values
are illustrative and must be benchmarked. Neither profile configures credentials
or institutional storage. No scheduler job has been submitted during setup.

`pilot/config.yaml` is the validated direct shared-server profile: local executor,
four cores and a conservative 16 GiB scheduling ceiling (mem_mb=16384). Start the
coordinator with nice -n 10 and monitor host/pilot memory and CPU. This profile
neither enforces an OS memory cap nor chooses a host. The completed pilot does
not authorize a new run merely because this profile is available.
