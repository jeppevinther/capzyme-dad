"""Optional sequential/genome-only orchestration configuration; external inputs."""
import csv
from pathlib import Path
import re

from configuration import checked_path, integer, keys
from genome_common import index_files
from reference_hashes import expected_digest, observed_file, check_manifest


def genome_settings(value, repo, output):
    if value is None:
        return None
    keys(value, ['aligner', 'index_prefix', 'reference', 'reference_sha256', 'decoys',
                 'no_decoys', 'verification_manifest', 'policy_files', 'threads', 'mem_mb', 'count_mem_mb',
                 'count_threads', 'mapq', 'exclude_duplicates', 'strand', 'seed'], 'genome')
    kind = value.get('aligner')
    if kind not in ('bowtie2', 'hisat2'):
        raise ValueError('genome.aligner must be bowtie2 or hisat2')
    prefix = str(checked_path(value.get('index_prefix'), 'genome.index_prefix', False))
    files = [str(p) for p in index_files(prefix, kind)]
    fasta = str(checked_path(value.get('reference'), 'genome.reference'))
    digest = expected_digest(value.get('reference_sha256'),'genome.reference_sha256')
    observed_file(fasta,digest)
    manifest=value.get('verification_manifest')
    if manifest is not None:
        manifest=str(checked_path(manifest,'genome.verification_manifest'))
        check_manifest(manifest,fasta,prefix,kind)
    decoys = value.get('decoys')
    decoys = str(checked_path(decoys, 'genome.decoys')) if decoys else None
    no_decoys = value.get('no_decoys', False)
    if type(no_decoys) is not bool or bool(decoys) == no_decoys:
        raise ValueError('Declare genome.decoys or genome.no_decoys: true explicitly')
    policies = value.get('policy_files', [])
    if not isinstance(policies, list):
        raise ValueError('genome.policy_files must be a list')
    policies = [str(checked_path(p, 'genome.policy_files')) for p in policies]
    strand = value.get('strand', 'both')
    if strand not in ('both', 'plus', 'minus'):
        raise ValueError('genome.strand must be plus, minus or both')
    duplicates = value.get('exclude_duplicates', False)
    if type(duplicates) is not bool:
        raise ValueError('genome.exclude_duplicates must be boolean')
    sources = [fasta, *files, *policies, *([manifest] if manifest else []), *([decoys] if decoys else [])]
    for source in sources:
        path = Path(source)
        if repo == path or repo in path.parents or output == path or output in path.parents:
            raise ValueError('Genome sources must be external and outside output_dir')
    return dict(aligner=kind, index_prefix=prefix, index_files=files, reference=fasta,
                reference_sha256=digest, verification_manifest=manifest, decoys=decoys, no_decoys=no_decoys,
                policy_files=policies, sources=sources, strand=strand, exclude_duplicates=duplicates,
                threads=integer(value.get('threads', 8), 'genome.threads', 1),
                mem_mb=integer(value.get('mem_mb', 32768), 'genome.mem_mb', 1),
                count_mem_mb=integer(value.get('count_mem_mb', 8192), 'genome.count_mem_mb', 1),
                count_threads=integer(value.get('count_threads', 2), 'genome.count_threads', 2),
                mapq=integer(value.get('mapq', 20), 'genome.mapq', 0, 255),
                seed=integer(value.get('seed', 0), 'genome.seed', 0, 2147483647))


def genome_only_settings(config, repo):
    keys(config, ['schema_version', 'samples', 'output_dir', 'genome', 'bigwigs'], 'genome-only config')
    if config.get('schema_version') != 1:
        raise ValueError('schema_version must be 1')
    out = checked_path(config.get('output_dir'), 'output_dir', False)
    if out == repo or repo in out.parents:
        raise ValueError('output_dir must be external')
    manifest = checked_path(config.get('samples'), 'samples')
    samples = {}
    with manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['sample_id', 'fastq']:
            raise ValueError('Manifest requires sample_id and fastq (already trimmed)')
        for row in reader:
            name = row['sample_id']
            if None in row or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', name) or name in samples or name in ('chromosome','strand','start','end'):
                raise ValueError('Invalid/duplicate/reserved sample name')
            fq = checked_path(row['fastq'], 'FASTQ')
            if out == fq or out in fq.parents or repo == fq or repo in fq.parents:
                raise ValueError('FASTQ must be external and outside output_dir')
            samples[name] = str(fq)
    if not samples or len(set(samples.values())) != len(samples):
        raise ValueError('Require unique samples and input files')
    if out == manifest or out in manifest.parents:
        raise ValueError('Manifest must be outside output_dir')
    genome = genome_settings(config.get('genome'), repo, out)
    if genome is None:
        raise ValueError('genome configuration is required')
    return dict(output=str(out), manifest=str(manifest), samples=samples, genome=genome)
