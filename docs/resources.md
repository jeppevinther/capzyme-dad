# Resources and shared-server execution

The [tutorial](getting_started.md) is deliberately tiny. Full experiments require
a separately reviewed resource/storage plan. First verify software, immutable
inputs, free disk space **and applicable quota**, available RAM and host load.

## Scheduling requests versus actual limits

`capzyme.py --cores N --memory-mb M` translates to Snakemake `--cores N
--resources mem_mb=M`. These are **total scheduling budgets** for concurrently
scheduled rules. Each rule also requests its own `threads` and `mem_mb`. A job
does not get a separate copy of the total budget: concurrent requests share it.
Choose a budget that can accommodate the largest job and the desired concurrency.
Snakemake may reduce a rule's threads to the available core budget.

The RNA streaming mapper runs Bowtie2, SAMtools filtering and sorting together;
use at least four total cores with its default four-thread job request. Sequential
mapping tools align, filter and sort in sequence within a job and reserve a helper
CPU slot. Genome threads/memory/count requests are set in the `genome` block;
existing-BAM requests are in `resources`. BAM preparations can run concurrently;
the cohort matrix is assembled after all BAMs are ready.

These settings are not operating-system enforcement. A tool can use more memory
than requested. A local executor does not enforce a hard total memory cap or CPU
affinity through these flags. Monitor the actual process tree and host reserve;
use site-approved cgroups/container/scheduler limits when hard enforcement is
needed. Do not treat an increased scheduling number as newly available physical
RAM. Keep numerical-library threads at one unless deliberately budgeted.

## Planning storage and memory

Use a reference- and dataset-specific reservation before full processing.
Human-scale HISAT2 index building can require substantial memory and scratch;
measure a representative approved run and retain sampled aggregate RSS separately
from maximum single-process RSS. Tiny fixtures establish correctness, not full-data
throughput or peak disk requirements.

Reserve output and temporary sorting space plus an explicit free-space floor.
Filesystem free space is not proof of project quota. Record any quota uncertainty
and obtain the storage owner's confirmation before a run that depends on it.
The four-core examples are scheduling examples, not universal production budgets.

## Disk use

Budget for source reads/references/indexes, compressed and uncompressed trimmed/
unmapped FASTQs, all-record and accepted BAMs, temporary SAM/SQLite/sort files,
unthresholded matrices, logs and retained attempts/history. Compressed raw FASTQ
bytes alone cannot predict those sizes. Generated-reference mode adds managed
FASTA and Bowtie2 indexes; genome mapping reuses an existing index, avoiding a new
large build. Multiple jobs multiply temporary space. Separate permanent references
from fresh per-analysis outputs and leave a free-space reserve.

## Shared servers

Use a site-approved compute host/allocation; ask whether direct local execution
is allowed. Do not assume that a login node permits heavy work. Host selection,
accounts, storage, credentials and any scheduler profile belong outside public
workflow code. The beginner launcher uses the local executor and is not a scheduler
submission client. Existing profile templates are not proof of operational SLURM
support at a particular site.

On an approved shared Linux server, lower CPU priority with `nice -n 10`:

```bash
nice -n 10 python "$CAPZYME_REPO/capzyme.py" run \
  --route sequential --config "$PRIVATE_CONFIG" --cores 32 --memory-mb 98304
```

Niceness affects CPU priority, **not** memory, disk usage or job isolation. Start
long work in an approved detached `screen`/`tmux` session so a terminal/network
interruption does not kill the session. For example, enter `screen -S capzyme`,
activate the environment and run the command inside it; `Ctrl-A`, then `D` detaches
and `screen -r capzyme` reconnects. This does not itself monitor or limit resources.

Use active process-tree RSS and host memory/CPU/load monitoring, not just a single
parent PID. Site-specific `free -h`, `uptime`, `df -h`, quota tools and process
monitors help inspect capacity; monitoring/stop thresholds belong to the external
execution plan. The retained shared-server policy reserves at least max(64 GiB,
10% of host RAM), stopping for sustained heavy CPU/load or breached memory/disk
reserves. That reserve targets large shared research servers, not the tiny laptop
tutorial. Stop and report pressure rather than moving work automatically. Record
commands, resolved versions, exit status, logs, actual resource use and locations.
