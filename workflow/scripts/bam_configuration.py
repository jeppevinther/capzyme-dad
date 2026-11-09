"""Configuration for the independent existing-BAM counting workflow."""
import csv
from pathlib import Path
import re

from bam_count import policy, sample_names
from configuration import checked_path, integer, keys
from core import sha256
from reference_hashes import expected_digest, observed_file


def load_bam_settings(config, repo):
    keys(config, ['schema_version', 'input_mode', 'samples', 'output_dir', 'reference',
                  'counting', 'comparison', 'resources', 'bigwigs'], 'config')
    if config.get('schema_version') != 1 or config.get('input_mode') != 'existing_bam':
        raise ValueError('Requires schema_version: 1 and input_mode: existing_bam')
    manifest = checked_path(config.get('samples'), 'samples')
    output = checked_path(config.get('output_dir'), 'output_dir', False)
    if output == repo or repo in output.parents:
        raise ValueError('output_dir must be outside the source checkout')
    ref = config.get('reference', {})
    keys(ref, ['fasta', 'fasta_sha256'], 'reference')
    fasta = checked_path(ref.get('fasta'), 'reference.fasta')
    digest = expected_digest(ref.get('fasta_sha256'), 'Reference SHA-256')
    observed_file(fasta,digest)
    samples = {}
    with manifest.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames != ['sample_id', 'bam']:
            raise ValueError('Manifest columns must be exactly sample_id and bam')
        for row in reader:
            name = row.get('sample_id')
            if None in row or name in samples:
                raise ValueError('Extra fields or duplicate sample ID')
            samples[name] = str(checked_path(row.get('bam'), 'BAM'))
    sample_names(list(samples), list(samples.values()))
    count = config.get('counting', {})
    keys(count, ['mapq', 'strand', 'exclude_duplicates', 'min_total_count'], 'counting')
    mapq, strand = count.get('mapq', 20), count.get('strand', 'both')
    duplicates, minimum = count.get('exclude_duplicates', False), count.get('min_total_count', 0)
    policy(mapq, duplicates, strand, minimum)
    comparison = config.get('comparison', {})
    keys(comparison, ['expected_matrix', 'expected_matrix_sha256', 'require_column_order',
                      'integer_tolerance', 'expected_features'], 'comparison')
    expected = None
    if comparison:
        expected = checked_path(comparison.get('expected_matrix'), 'comparison.expected_matrix')
        if sha256(expected) != comparison.get('expected_matrix_sha256'):
            raise ValueError('Historical matrix SHA-256 mismatch')
        if comparison.get('require_column_order', True) is not True or comparison.get('integer_tolerance', 0) != 0:
            raise ValueError('Comparison requires exact column order and zero integer tolerance')
        if 'expected_features' in comparison:
            integer(comparison['expected_features'], 'expected_features', 0)
    resource = config.get('resources', {})
    keys(resource, ['cores', 'mem_mb', 'filter_threads', 'filter_mem_mb', 'count_threads', 'count_mem_mb'], 'resources')
    for key, default in [('cores', 4), ('mem_mb', 16384), ('filter_threads', 1),
                         ('filter_mem_mb', 1024), ('count_threads', 2), ('count_mem_mb', 8192)]:
        integer(resource.get(key, default), 'resources.'+key, 1)
    if resource.get('filter_threads', 1) != 1:
        raise ValueError('BAM preparation uses one thread per job')
    inputs = [manifest, fasta, *(Path(b) for b in samples.values()), *([expected] if expected else [])]
    for path in inputs:
        if path == output or output in path.parents:
            raise ValueError('Inputs cannot reside inside output_dir')
        if path != manifest and (path == repo or repo in path.parents):
            raise ValueError('Research inputs must remain outside checkout')
    return dict(output=str(output), manifest=str(manifest), samples=samples, fasta=str(fasta),
                reference_sha256=digest, mapq=mapq, strand=strand, exclude_duplicates=duplicates,
                minimum=minimum, comparison=dict(comparison), resources=dict(resource))
