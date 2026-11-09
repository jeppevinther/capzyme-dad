"""Milestone 1 configuration contract; paths resolve against the working directory."""
import csv
from pathlib import Path
import re
from fixed_reference import check_fixed_hashes, fixed_inputs


def integer(value, label, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f'{label} must be an integer in [{minimum}, {maximum}]')
    return value


def keys(mapping, allowed, label):
    if not isinstance(mapping, dict):
        raise ValueError(f'{label} must be a mapping')
    extra = set(mapping) - set(allowed)
    if extra:
        raise ValueError(f'Unknown {label} fields: {sorted(extra)}')


def checked_path(value, label, required_file=True):
    if not isinstance(value, str) or not value or any(c in value for c in '\n\r\t{}'):
        raise ValueError(f'{label} requires a nonempty path without control characters or braces')
    path = Path(value).expanduser().resolve()
    if required_file and not path.is_file():
        raise ValueError(f'{label} does not exist: {path}')
    return path


def load_settings(config, repo):
    keys(config, ['schema_version', 'samples', 'output_dir', 'reference', 'preprocessing', 'mapping', 'counting', 'genome', 'bigwigs'], 'config')
    if config.get('schema_version') != 1:
        raise ValueError('schema_version must be 1; provide --configfile with a milestone 1 configuration')
    output = checked_path(config.get('output_dir'), 'output_dir', False)
    if output == repo or repo in output.parents:
        raise ValueError('output_dir must be outside the source checkout')
    manifest = checked_path(config.get('samples'), 'samples')
    ref = config.get('reference', {})
    mode = ref.get('mode', 'generated')
    if mode not in ('generated', 'fixed'):
        raise ValueError('reference.mode must be generated or fixed')
    fixed = None
    if mode == 'fixed':
        keys(ref, ['mode', 'fasta', 'fasta_sha256', 'index_prefix', 'index_format', 'index_sha256'], 'reference')
        fixed = dict(ref)
        fixed['fasta'] = str(checked_path(ref.get('fasta'), 'reference.fasta'))
        fixed['index_prefix'] = str(checked_path(ref.get('index_prefix'), 'reference.index_prefix', False))
        check_fixed_hashes(fixed)  # Reject missing/changed artifacts before scheduling reads.
        _, index_files = fixed_inputs(fixed)
        reference_sources = [fixed['fasta'], *index_files]
        ncrna, viral, k, drops = None, [], None, []
    else:
        keys(ref, ['mode', 'ncrna_fasta', 'viral_fastas', 'prefix_length', 'drop_biotypes'], 'reference')
        ncrna = checked_path(ref.get('ncrna_fasta'), 'reference.ncrna_fasta')
        viral = ref.get('viral_fastas', [])
        if not isinstance(viral, list):
            raise ValueError('reference.viral_fastas must be a list')
        viral = [checked_path(x, 'viral_fastas') for x in viral]
        k = integer(ref.get('prefix_length', 30), 'prefix_length', 1)
        drops = ref.get('drop_biotypes', [])
        if not isinstance(drops, list) or any(not isinstance(x, str) or not re.fullmatch(r'[A-Za-z0-9_+-]+', x) for x in drops):
            raise ValueError('drop_biotypes must be a list of biotype names')
        reference_sources = [str(ncrna), *(str(x) for x in viral)]
    pre = config.get('preprocessing', {})
    keys(pre, ['adapter', 'quality', 'minimum_length', 'nextseq_trim'], 'preprocessing')
    if not isinstance(pre.get('adapter'), str) or not re.fullmatch('[ACGTNacgtn]+', pre['adapter']):
        raise ValueError('preprocessing.adapter must be explicitly supplied as a DNA sequence')
    for field, lower, upper in [('quality', 0, 93), ('minimum_length', 1, None), ('nextseq_trim', 0, 93)]:
        if field=='nextseq_trim' and field in pre and pre[field] is None: continue
        integer(pre.get(field), f'preprocessing.{field}', lower, upper)
    mapping = config.get('mapping', {})
    keys(mapping, ['norc', 'preset', 'seed', 'threads', 'keep_diagnostic_sam'], 'mapping')
    norc = mapping.get('norc', True)
    if type(norc) is not bool:
        raise ValueError('mapping.norc must be a boolean')
    preset = mapping.get('preset', 'very-sensitive')
    if preset != 'very-sensitive':
        raise ValueError('Only the established very-sensitive preset is supported in milestone 1')
    diagnostic = mapping.get('keep_diagnostic_sam', False)
    if type(diagnostic) is not bool:
        raise ValueError('mapping.keep_diagnostic_sam must be a boolean')
    seed = integer(mapping.get('seed', 0), 'mapping.seed', 0, 2147483647)
    threads = integer(mapping.get('threads', 4), 'mapping.threads', 1)
    counting = config.get('counting', {})
    keys(counting, ['mapq', 'strand', 'exclude_duplicates', 'min_total_count'], 'counting')
    mapq = integer(counting.get('mapq', 20), 'counting.mapq', 0, 255)
    minimum = integer(counting.get('min_total_count', 0), 'counting.min_total_count', 0)
    strand = counting.get('strand', 'both')
    if strand not in ('plus', 'minus', 'both'):
        raise ValueError('counting.strand must be plus, minus or both')
    exclude_duplicates = counting.get('exclude_duplicates', False)
    if type(exclude_duplicates) is not bool:
        raise ValueError('counting.exclude_duplicates must be a boolean')
    samples = {}
    with manifest.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames != ['sample_id', 'fastq']:
            raise ValueError('Manifest columns must be exactly sample_id and fastq (single-end)')
        for row in reader:
            name = row.get('sample_id', '')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name) or name in samples or name in ('chromosome', 'strand', 'start', 'end'):
                raise ValueError(f'Invalid, duplicate, or reserved sample ID: {name!r}')
            if None in row:
                raise ValueError(f'Extra manifest fields for {name}')
            samples[name] = str(checked_path(row.get('fastq'), f'FASTQ for {name}'))
    if not samples:
        raise ValueError('Manifest must contain at least one sample')
    if len(set(samples.values())) != len(samples):
        raise ValueError('The same FASTQ must not be assigned to multiple samples')
    for path in [manifest, *(Path(x) for x in reference_sources), *(Path(x) for x in samples.values())]:
        if path == output or output in path.parents:
            raise ValueError('Inputs cannot reside inside output_dir (prevents input/output collisions)')
        if path != manifest and (path == repo or repo in path.parents):
            raise ValueError('Sequencing/reference files must remain outside the source checkout')
    from genome_configuration import genome_settings
    genome = genome_settings(config.get('genome'), repo, output)
    return dict(output=str(output), manifest=str(manifest), samples=samples, ncrna=(str(ncrna) if ncrna else None),
                viral=[str(x) for x in viral], k=k, drops=drops, pre=pre,
                reference_mode=mode, fixed=fixed, reference_sources=reference_sources,
                norc=norc, preset=preset, seed=seed, threads=threads, diagnostic=diagnostic,
                mapq=mapq, minimum=minimum, strand=strand,
                flags=2308 | (1024 if exclude_duplicates else 0), genome=genome)
