"""Strict version-2 standalone settings adapter. Scientific engines are unchanged."""
import argparse
import csv
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

from core import collapse, sha256
from genome_common import dictionary, executable, index_files, inspector_tool

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'workflow/scripts'


def keys(value, allowed, label):
    if not isinstance(value, dict) or set(value) - set(allowed):
        raise ValueError(f'Unsupported fields in {label}: {set(value) - set(allowed) if isinstance(value, dict) else value}')


def load(path):
    if not path:
        return {}, Path.cwd()
    path = Path(path).expanduser().resolve()
    # JSON is supported without PyYAML; YAML requires only PyYAML, never Snakemake.
    text = path.read_text()
    def unique_pairs(pairs):
        result = {}
        for key, item in pairs:
            if key in result: raise ValueError('Duplicate settings key: ' + str(key))
            result[key] = item
        return result
    try:
        value = json.loads(text, object_pairs_hook=unique_pairs)
    except json.JSONDecodeError:
        import yaml
        class UniqueLoader(yaml.SafeLoader):
            pass
        def unique(loader, node):
            result = {}
            for k, v in node.value:
                key = loader.construct_object(k)
                if key in result:
                    raise ValueError('Duplicate settings key: ' + str(key))
                result[key] = loader.construct_object(v)
            return result
        UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique)
        try: value = yaml.load(text, Loader=UniqueLoader)
        except yaml.YAMLError as exc: raise ValueError('Invalid YAML: ' + str(exc)) from exc
    keys(value, ['schema_version', 'experiment_name', 'samples', 'output_dir', 'input',
                 'references', 'trimming', 'mapping', 'counting', 'bigwig'], 'settings')
    if type(value.get('schema_version')) is not int or value['schema_version'] != 2:
        raise ValueError('toolbox settings require schema_version: 2; use capzyme.py for version 1')
    for name in ('experiment_name', 'samples', 'output_dir'):
        if name in value and not isinstance(value[name], str):
            raise ValueError(name + ' must be a string')
    keys(value.get('input', {}), ['type', 'matrix', 'assert_unthresholded'], 'input')
    if 'type' in value.get('input', {}) and value['input']['type'] not in ('raw_fastq', 'trimmed_fastq', 'bam', 'counts'):
        raise ValueError('Unsupported input.type')
    keys(value.get('references', {}), ['rna', 'genome'], 'references')
    for space, ref in value.get('references', {}).items():
        keys(ref, ['fasta', 'chrom_sizes', 'index', 'decoys_file', 'no_decoys'], 'references.' + space)
        keys(ref.get('index', {}), ['action', 'prefix', 'verification_manifest'], space + '.index')
        if ref.get('index', {}).get('action', 'use_existing') != 'use_existing':
            raise ValueError('Only use_existing is supported here; build_new is an explicit native preparation operation')
    keys(value.get('mapping', {}), ['rna', 'genome'], 'mapping')
    for space, mapping in value.get('mapping', {}).items():
        allowed = ['aligner', 'preset', 'seed', 'minimum_mapq', 'exclude_duplicates', 'threads']
        allowed += ['allow_reverse_complement'] if space == 'rna' else ['uncertain_cigar']
        keys(mapping, allowed, 'mapping.' + space)
        if mapping.get('preset', 'very-sensitive') != 'very-sensitive':
            raise ValueError('Only the validated very-sensitive preset is supported')
        if mapping.get('aligner', 'bowtie2' if space == 'rna' else 'hisat2') not in (('bowtie2',) if space == 'rna' else ('hisat2', 'bowtie2')):
            raise ValueError('Unsupported aligner for ' + space)
    keys(value.get('trimming', {}), ['adapter', 'quality', 'minimum_length', 'nextseq_trim'], 'trimming')
    keys(value.get('counting', {}), ['rna_strands', 'genome_strands', 'minimum_total_count', 'minimum_mapq', 'exclude_duplicates'], 'counting')
    keys(value.get('bigwig', {}), ['rna', 'genome'], 'bigwig')
    for space, bw in value.get('bigwig', {}).items():
        keys(bw, ['strands', 'normalization', 'cohort', 'emit_samples', 'depths', 'matrix_complete_both_strands',
                  'deseq2_type', 'controls', 'min_control_features', 'size_factors', 'factor_kind', 'factor_metadata'], 'bigwig.' + space)
    # Validate every supported scalar, even in sections not consumed by this operation.
    def validate(node):
        for k, v in node.items():
            if isinstance(v, dict):
                validate(v)
            elif v is None and k in ('nextseq_trim','verification_manifest'):
                continue
            elif k in ('allow_reverse_complement', 'exclude_duplicates', 'no_decoys', 'assert_unthresholded', 'matrix_complete_both_strands'):
                if type(v) is not bool: raise ValueError(k + ' must be boolean')
            elif k in ('seed', 'minimum_mapq', 'threads', 'minimum_total_count', 'quality', 'minimum_length', 'nextseq_trim', 'min_control_features'):
                if type(v) is not int or v < 0: raise ValueError(k + ' must be a nonnegative integer')
                if k == 'minimum_mapq' and v > 255: raise ValueError('minimum_mapq must be <=255')
                if k in ('threads', 'minimum_length') and v < 1: raise ValueError(k + ' must be positive')
                if k == 'min_control_features' and v < 2: raise ValueError('min_control_features must be >=2')
            elif k in ('cohort', 'emit_samples'):
                if not isinstance(v, list) or not v or any(not isinstance(x, str) or not x for x in v) or len(set(v)) != len(v):
                    raise ValueError(k + ' requires unique nonempty string sample IDs')
            elif not isinstance(v, str):
                raise ValueError(k + ' must be a string')
            elif k in ('strands', 'rna_strands', 'genome_strands') and v not in ('plus', 'minus', 'both'):
                raise ValueError('Invalid strand choice')
            elif k == 'normalization' and v not in ('raw', 'rpm', 'deseq2', 'controls', 'supplied'):
                raise ValueError('Unsupported normalization')
            elif k == 'deseq2_type' and v not in ('ratio', 'poscounts'):
                raise ValueError('Unsupported DESeq2 estimator')
            elif k == 'factor_kind' and v not in ('size-factor', 'multiplier'):
                raise ValueError('Unsupported factor_kind')
            elif k == 'uncertain_cigar' and v not in ('reject', 'error'):
                raise ValueError('Unsupported CIGAR policy')
    validate({k: v for k, v in value.items() if k != 'schema_version'})
    return value, path.parent


def absolute(value, base):
    p = Path(value).expanduser()
    return str((base / p).resolve()) if not p.is_absolute() else str(p.resolve())


def samples(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        allowed = {'sample_id', 'input_file', 'experiment', 'treatment', 'biological_source', 'batch', 'display_label'}
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or set(reader.fieldnames) - allowed or 'sample_id' not in reader.fieldnames:
            raise ValueError('Sample sheet needs sample_id and supported metadata/input_file columns')
        rows = list(reader)
    ids = [r['sample_id'] for r in rows]
    if not rows or len(set(ids)) != len(ids) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', x) or x in ('start', 'end', 'chromosome', 'strand') for x in ids):
        raise ValueError('Invalid/duplicate sample IDs')
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('Malformed sample row')
    for row in rows:
        if row.get('input_file'):
            row['input_file'] = absolute(row['input_file'], Path(path).parent)
    return rows


from reference_hashes import check_manifest


def parser():
    p = argparse.ArgumentParser(prog='capzyme', description=__doc__)
    sub = p.add_subparsers(dest='operation', required=True)
    for op in ('trim', 'rna-align', 'genome-align', 'count', 'bigwig'):
        q = sub.add_parser(op, help='Standalone ' + op + '; --preview prints resolved settings and commands')
        q.add_argument('--settings', help='Strict version-2 YAML/JSON; relative paths bind to this file')
        q.add_argument('--preview', action='store_true')
        q.add_argument('--output-dir', help='New external output directory')
        q.add_argument('--samples', help='TSV: sample_id, input_file plus optional metadata')
        q.add_argument('--sample', help='One sample for trim/alignment')
        q.add_argument('--input-file', help='Explicit FASTQ; overrides selected sample input')
        q.add_argument('--space', choices=['rna', 'genome'], default=None)
        q.add_argument('--reference'); q.add_argument('--chrom-sizes')
        q.add_argument('--matrix'); q.add_argument('--bam', action='append')
        q.add_argument('--sample-name', action='append')
        q.add_argument('--bam-manifest')
        q.add_argument('--index'); q.add_argument('--verification-manifest')
        q.add_argument('--aligner', choices=['bowtie2', 'hisat2'])
        q.add_argument('--threads', type=int); q.add_argument('--seed', type=int)
        q.add_argument('--mapq', type=int)
        for name in ('exclude-duplicates', 'allow-reverse-complement', 'assert-unthresholded', 'matrix-complete-both-strands'):
            q.add_argument('--' + name, action=argparse.BooleanOptionalAction, default=None)
        q.add_argument('--no-decoys', action='store_true', default=None, help='Explicitly declare a reference without decoys')
        q.add_argument('--decoys'); q.add_argument('--uncertain-cigar', choices=['reject', 'error'])
        q.add_argument('--strand', choices=['plus', 'minus', 'both'])
        q.add_argument('--min-total-count', type=int)
        q.add_argument('--adapter'); q.add_argument('--quality', type=int)
        q.add_argument('--minimum-length', type=int)
        group=q.add_mutually_exclusive_group()
        group.add_argument('--nextseq-trim', type=int, help='Two-color trimming threshold; zero still enables treatment')
        group.add_argument('--no-nextseq-trim', dest='nextseq_trim', action='store_const', const=False, help='Omit NextSeq-specific treatment; ordinary quality trimming remains active')
        q.add_argument('--normalization', choices=['raw', 'rpm', 'deseq2', 'controls', 'supplied'])
        q.add_argument('--cohort', nargs='+'); q.add_argument('--emit-samples', nargs='+')
        q.add_argument('--depths'); q.add_argument('--deseq2-type', choices=['ratio', 'poscounts'])
        q.add_argument('--controls'); q.add_argument('--min-control-features', type=int)
        q.add_argument('--size-factors'); q.add_argument('--factor-kind', choices=['size-factor', 'multiplier'])
        q.add_argument('--factor-metadata'); q.add_argument('--rscript')
    q = sub.add_parser('collapse', help='Unchanged legacy prefix-collapse implementation')
    q.add_argument('--reference', required=True); q.add_argument('--output-dir', required=True)
    q.add_argument('--prefix-length', type=int, default=30); q.add_argument('--drop-biotype', action='append', default=[])
    q = sub.add_parser('index-record', help='Record existing native index hashes/dictionary; newly observed unless trusted manifest supplied')
    q.add_argument('--reference', required=True); q.add_argument('--index', required=True)
    q.add_argument('--aligner', required=True, choices=['bowtie2', 'hisat2'])
    q.add_argument('--trusted-manifest'); q.add_argument('--output-dir', required=True)
    return p


def resolve(a):
    cfg, base = load(a.settings)
    op = a.operation
    space = a.space or ('genome' if op == 'genome-align' else 'rna')
    if op == 'rna-align' and space != 'rna' or op == 'genome-align' and space != 'genome':
        raise ValueError('Alignment operation and --space disagree')
    ref = cfg.get('references', {}).get(space, {})
    mapping = cfg.get('mapping', {}).get(space, {})
    bw = cfg.get('bigwig', {}).get(space, {})
    count = cfg.get('counting', {})
    trim = cfg.get('trimming', {})
    used = set()
    def get(name, value=None, default=None, path=False):
        used.add(name)
        cli = getattr(a, name, None)
        v = cli if cli is not None else (value if value is not None else default)
        if v is not None and path:
            v = absolute(v, Path.cwd() if cli is not None else base)
        return v
    out = get('output_dir', cfg.get('output_dir'), path=True)
    if not out: raise ValueError('Supply --output-dir or output_dir')
    command = []
    rows = []
    sheet = get('samples', cfg.get('samples'), path=True)
    if sheet: rows = samples(sheet)
    def required(v, label):
        if v is None or v == []: raise ValueError('Required: ' + label)
        return v
    def option(flag, value):
        if value is not None:
            command.append('--' + flag)
            command.extend(map(str, value if isinstance(value, list) else [value]))
    def truth(flag, value):
        if value: command.append('--' + flag)
    policies = {'excluded_flags': 2308, 'QC_failure_records': 'retained', 'layout': 'single-end'}
    if op in ('trim', 'rna-align', 'genome-align'):
        selected = get('sample')
        input_file = get('input_file', path=True)
        if not input_file:
            row = next((r for r in rows if r['sample_id'] == selected), None)
            if not row or not row.get('input_file'): raise ValueError('Select --sample from sheet or supply --input-file')
            input_file = row['input_file']
        typ = cfg.get('input', {}).get('type')
        if typ and typ != ('raw_fastq' if op == 'trim' else 'trimmed_fastq'):
            raise ValueError('input.type must explicitly match this operation (raw_fastq for trim; trimmed_fastq for alignment)')
    if op == 'trim':
        command = ['cutadapt', '-a', required(get('adapter', trim.get('adapter')), 'adapter'),
                   '-q', str(get('quality', trim.get('quality'), 20)),
                   '-m', str(get('minimum_length', trim.get('minimum_length'), 15)),
                   '--cores', str(get('threads', default=1)), '--json', str(Path(out)/'cutadapt.json'),
                   '-o', str(Path(out)/'trimmed.fastq.gz'), input_file]
        used.add('nextseq_trim')
        nextseq = a.nextseq_trim if a.nextseq_trim is not None else trim.get('nextseq_trim',20)
        if nextseq is not None and nextseq is not False:
            if type(nextseq) is not int or not 0<=nextseq<=93: raise ValueError('nextseq_trim must be null or an integer from 0 to 93')
            command[1:1]=['--nextseq-trim',str(nextseq)]
    else:
        reference = get('reference', ref.get('fasta'), path=True)
        if op in ('rna-align', 'genome-align'):
            kind = get('aligner', mapping.get('aligner'), 'bowtie2' if space == 'rna' else 'hisat2')
            if op == 'rna-align' and kind != 'bowtie2': raise ValueError('RNA alignment requires Bowtie2')
            index = required(get('index', ref.get('index', {}).get('prefix'), path=True), 'existing index')
            manifest = get('verification_manifest', ref.get('index', {}).get('verification_manifest'), path=True)
            required(reference, 'finished reference FASTA')
            if manifest: check_manifest(manifest, reference, index, kind)
            engine = 'rna_handoff.py' if op == 'rna-align' else kind + '_genome.py'
            command = [sys.executable, str(SCRIPTS/engine)] + ([] if op == 'rna-align' else ['map'])
            option('fastq', input_file); option('reference', reference); option('index', index)
            option('threads', get('threads', mapping.get('threads'), 1))
            option('seed', get('seed', mapping.get('seed'), 0))
            option('mapq', get('mapq', mapping.get('minimum_mapq'), 20))
            truth('exclude-duplicates', get('exclude_duplicates', mapping.get('exclude_duplicates'), False))
            if op == 'rna-align':
                truth('norc', not get('allow_reverse_complement', mapping.get('allow_reverse_complement'), False))
                policies['orientation'] = 'RNA reverse complement explicitly resolved via --norc presence'
            else:
                decoys = get('decoys', ref.get('decoys_file'), path=True)
                no_decoys = get('no_decoys', ref.get('no_decoys'), False)
                if bool(decoys) == bool(no_decoys): raise ValueError('Declare exactly one of decoys_file / no_decoys')
                option('decoys', decoys); truth('no-decoys', no_decoys)
                option('uncertain-cigar', get('uncertain_cigar', mapping.get('uncertain_cigar'), 'reject'))
                policies.update(preset='very-sensitive', orientations='both', soft_clipping='disabled',
                                CIGAR='terminal M/= /X; no S/H/P; internal I/D/N allowed', sort_memory='256M per worker')
        else:
            command = [sys.executable, str(SCRIPTS/('count_5ends.py' if op == 'count' else 'bigwig2.py'))]
            if op == 'count':
                required(reference, 'reference FASTA')
                option('reference', reference)
                strand = get('strand', count.get(space + '_strands'), 'both')
                option('strand', strand)
                option('mapq', get('mapq', count.get('minimum_mapq'), 20))
                truth('exclude-duplicates', get('exclude_duplicates', count.get('exclude_duplicates'), False))
                option('min-total-count', get('min_total_count', count.get('minimum_total_count'), 0))
                bams = get('bam') or [r.get('input_file') for r in rows]
                names = get('sample_name') or [r['sample_id'] for r in rows]
                required(bams, 'BAM inputs'); required(names, 'ordered sample names')
                if any(not b for b in bams) or len(bams) != len(names): raise ValueError('One sample name required per BAM')
                if cfg.get('input', {}).get('type', 'bam') != 'bam': raise ValueError('count settings require input.type: bam')
                for name in names: option('sample-name', name)
                option('output', str(Path(out)/'counts.tsv'))
                command += [absolute(b, Path.cwd()) for b in bams]
                policies['coordinates'] = '1-based start=end; reference-relative strand'
            else:
                matrix = get('matrix', cfg.get('input', {}).get('matrix'), path=True)
                bams = get('bam')
                manifest = get('bam_manifest', path=True)
                if sum(bool(x) for x in (matrix, bams, manifest)) != 1:
                    raise ValueError('Supply exactly one matrix, repeated --bam, or --bam-manifest')
                typ = cfg.get('input', {}).get('type')
                if typ and typ != ('counts' if matrix else 'bam'):
                    raise ValueError('bigwig input.type must match counts or bam input')
                option('matrix', matrix)
                if matrix:
                    if not get('assert_unthresholded', cfg.get('input', {}).get('assert_unthresholded'), False):
                        raise ValueError('Explicit --assert-unthresholded or input.assert_unthresholded: true required')
                    truth('assert-unthresholded', True)
                else:
                    for bam in bams or []: option('bam', absolute(bam, Path.cwd()))
                    option('bam-manifest', manifest)
                    for name in get('sample_name') or []: option('sample-name', name)
                    option('mapq', get('mapq', count.get('minimum_mapq'), 20))
                    truth('exclude-duplicates', get('exclude_duplicates', count.get('exclude_duplicates'), False))
                sizes = get('chrom_sizes', ref.get('chrom_sizes'), path=True)
                if bool(reference) == bool(sizes): raise ValueError('Supply exactly one FASTA/chrom_sizes')
                option('reference', reference); option('chrom-sizes', sizes)
                option('normalization', required(get('normalization', bw.get('normalization')), 'normalization'))
                option('cohort', get('cohort', bw.get('cohort'), [r['sample_id'] for r in rows] or None))
                option('emit-samples', get('emit_samples', bw.get('emit_samples')))
                option('strands', get('strand', bw.get('strands'), 'both'))
                for name in ('depths', 'controls', 'size_factors', 'factor_metadata'):
                    option(name.replace('_', '-'), get(name, bw.get(name), path=True))
                for name in ('deseq2_type', 'factor_kind', 'min_control_features'):
                    option(name.replace('_', '-'), get(name, bw.get(name)))
                option('rscript', get('rscript'))
                truth('matrix-complete-both-strands', get('matrix_complete_both_strands', bw.get('matrix_complete_both_strands'), False))
                # Use the existing parser for option combinations; execute() owns scientific validation.
                policies['coordinates'] = '1-based p -> BigWig [p-1,p); nonnegative both strands'
        if op != 'count': option('output-dir', out)
    common = {'operation', 'settings', 'preview', 'space'}
    for name, value in vars(a).items():
        if name not in used | common and value is not None:
            raise ValueError('--' + name.replace('_', '-') + ' is unsupported for ' + op)
    # Validate effective numeric CLI overrides too, before a tool creates outputs.
    for flag in ('--mapq', '--seed', '--threads', '--cores', '--min-total-count', '-q', '-m', '--nextseq-trim'):
        if flag in command:
            number = int(command[command.index(flag)+1])
            if number < (1 if flag in ('--threads', '--cores', '-m') else 0) or flag == '--mapq' and number > 255:
                raise ValueError('Invalid effective ' + flag)
    if '--exclude-duplicates' in command: policies['excluded_flags'] |= 1024
    engine_parameters = None
    if op == 'bigwig':
        from bigwig2 import parser as bigwig_parser, validate_options
        engine_args = bigwig_parser().parse_args(command[2:])
        validate_options(engine_args)
        engine_parameters = vars(engine_args)
        if engine_parameters['min_control_features'] is None: engine_parameters['min_control_features']=10
        # Reject invalid option combinations at preview time without touching data.
        x = engine_parameters
        if x['matrix'] and x['normalization'] == 'rpm' and not x['depths'] and not x['matrix_complete_both_strands']: raise ValueError('RPM needs explicit depths or complete-both-strands assertion')
    return {'operation': op, 'space': space, 'command': command, 'command_shell': shlex.join(command),
            'resolved_engine_parameters': engine_parameters,
            'output_dir': out, 'policies': policies, 'settings': str(Path(a.settings).resolve()) if a.settings else None,
            'settings_sha256': sha256(a.settings) if a.settings else None,
            'samples': rows, 'precedence': 'explicit CLI > settings > documented adapter defaults'}


def run(a):
    if a.operation in ('collapse', 'index-record'):
        out = Path(a.output_dir).expanduser().resolve()
        if out.exists() or out == ROOT or ROOT in out.parents: raise ValueError('Choose a new external output directory')
        source = Path(a.reference).expanduser().resolve()
        before = sha256(source)
        if a.operation == 'collapse' and a.prefix_length < 1: raise ValueError('prefix-length must be positive')
        if a.operation == 'index-record' and a.trusted_manifest:
            check_manifest(a.trusted_manifest, source, a.index, a.aligner)
        out.mkdir(parents=True)
        report = {'status': 'running', 'operation': a.operation, 'parameters': vars(a), 'reference_sha256': before}
        try:
            if a.operation == 'collapse':
                report['collapse'] = collapse(source, out/'collapsed', a.prefix_length, a.drop_biotype, ROOT/'legacy_scripts/merge_by_prefix_v2.py')
            else:
                parts = index_files(a.index, a.aligner)
                pins = {str(p)[len(str(Path(a.index).resolve())):]: sha256(p) for p in parts}
                with (out/'commands.jsonl').open('w') as log:
                    sizes = dictionary(a.index, a.aligner, out, log, source)
                if pins != {str(p)[len(str(Path(a.index).resolve())):]: sha256(p) for p in parts}: raise ValueError('Index changed during recording')
                manifest = {'schema_version': 1, 'aligner': a.aligner, 'reference_sha256': before, 'parts': pins,
                            'dictionary': list(sizes.items()), 'verification': 'FASTA/index identifiers and lengths; sequence identity requires trustworthy pins/build provenance',
                            'trust': 'matched supplied trusted manifest' if a.trusted_manifest else 'newly observed hashes; not pre-existing trusted evidence'}
                (out/'reference_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
                report['software'] = executable(inspector_tool(a.index, a.aligner))
            if sha256(source) != before: raise ValueError('Reference changed')
            report['status'] = 'complete'
        except Exception as exc:
            report.update(status='failed', error=str(exc)); raise
        finally:
            (out/'methods.json').write_text(json.dumps(report, indent=2)+'\n')
        return
    plan = resolve(a)
    print(json.dumps(plan, indent=2), flush=True)
    if a.preview: return
    out = Path(plan['output_dir'])
    if out.exists() or out == ROOT or ROOT in out.parents: raise ValueError('Choose a new external output directory')
    # Mapping/BigWig own their output creation. Trimming/counting need a parent.
    if a.operation in ('trim', 'count'): out.mkdir(parents=True)
    started = time.time()
    try:
        if a.operation == 'trim':
            plan['software'] = executable('cutadapt')
            with (out/'cutadapt.log').open('w') as log:
                result = subprocess.run(plan['command'], stdout=log, stderr=subprocess.STDOUT, check=False)
        else:
            result = subprocess.run(plan['command'], check=False)
        plan.update(exit_status=result.returncode, elapsed_seconds=time.time()-started)
        if result.returncode: raise subprocess.CalledProcessError(result.returncode, plan['command'])
    finally:
        if out.is_dir():
            plan.setdefault('exit_status', 1)
            plan['adapter_sha256'] = sha256(__file__)
            (out/'methods.json').write_text(json.dumps(plan, indent=2)+'\n')


def main():
    p = parser()
    try: run(p.parse_args())
    except (ValueError, OSError, ImportError, subprocess.SubprocessError) as exc:
        p.exit(1, 'toolbox: ' + str(exc) + '\n')
