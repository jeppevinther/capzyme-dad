"""Workflow I/O coordination and exclusive read accounting; no counting engine."""
from collections import Counter
from contextlib import closing
import json
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
import tempfile
import datetime
import uuid

from bam_count import command
from core import sha256
from genome_common import fastq_records
from rna_handoff import classify


def publish_stage(root, space, sample, args, log):
    """Run a standalone tool in a fresh attempt and publish only after success.

    Failed attempts and previous workflow-owned results are preserved. No
    historical input or unmanaged directory is removed to satisfy CLI freshness.
    """
    root = Path(root).resolve()
    parent = root/'.attempts'/space/sample
    parent.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=parent))
    result = attempt/'result'
    command([*args, '--output-dir', str(result)], log)
    meta = json.loads((result/'provenance.json').read_text())
    if meta['status'] != 'complete' or not meta['inputs_unchanged']:
        raise ValueError('Standalone stage did not complete with unchanged inputs')
    destination = root/space/sample
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        raise ValueError('Refuse workflow output symlink')
    if destination.exists():
        marker = destination/'.workflow-owned'
        # Snakemake pre-creates output parents. An empty directory contains no
        # user result; preserve it in history just as an owned previous result.
        if any(destination.iterdir()) and (not marker.is_file() or marker.read_text() != str(destination)):
            raise ValueError('Refuse replacing unmanaged stage directory')
        history = root/'.history'/space/sample
        history.mkdir(parents=True, exist_ok=True)
        destination.rename(history/str(uuid.uuid4()))
    (result/'.workflow-owned').write_text(str(destination))
    (result/'publication.json').write_text(json.dumps({
        'standalone_output_dir': str(result), 'published_output_dir': str(destination),
        'standalone_provenance_sha256': sha256(result/'provenance.json'),
        'published_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'policy': 'Successful fresh standalone attempt; original tool provenance preserved'
    }, indent=2)+'\n')
    result.rename(destination)
    return destination


def bam_records(path, log):
    args = ['samtools', 'view', str(path)]
    log.write(json.dumps({'command': args})+'\n'); log.flush()
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log, text=True)
    try:
        for line in proc.stdout:
            yield line.rstrip('\n').split('\t')
        code = proc.wait()
        if code:
            raise subprocess.CalledProcessError(code, args)
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate(); proc.wait()


def _json_object(path):
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f'Missing or invalid integrity record: {path}') from exc
    if not isinstance(value, dict):
        raise ValueError(f'Integrity record must be an object: {path}')
    return value


def _digest(value, label):
    if not isinstance(value, str) or not re.fullmatch('[a-f0-9]{64}', value):
        raise ValueError(f'Missing or invalid SHA-256: {label}')
    return value


def verify_stage(directory, handoff=False):
    """Bind consumed artifacts to complete provenance, including publication.

    Output hashes in standalone provenance are relative to its original output
    directory. A publication receipt proves their relocation without rewriting
    the original provenance or searching by basename among unrelated inputs.
    """
    directory = Path(directory).resolve()
    provenance = directory/'provenance.json'
    meta = _json_object(provenance)
    if meta.get('status') != 'complete' or meta.get('inputs_unchanged') is not True:
        raise ValueError(f'Stage provenance is not complete/input-verified: {directory}')
    parameters = meta.get('parameters')
    if (not isinstance(parameters, dict) or not isinstance(parameters.get('output_dir'), str)
            or not parameters['output_dir']):
        raise ValueError(f'Missing stage output-directory provenance: {directory}')
    original = Path(parameters['output_dir'])
    original = original.resolve()
    pins = {str(provenance): sha256(provenance)}
    publication = directory/'publication.json'
    if publication.exists() or original != directory:
        receipt = _json_object(publication)
        if (receipt.get('standalone_output_dir') != str(original) or
                receipt.get('published_output_dir') != str(directory) or
                receipt.get('standalone_provenance_sha256') != pins[str(provenance)]):
            raise ValueError(f'Publication/provenance connection mismatch: {directory}')
        pins[str(publication)] = sha256(publication)
    outputs = meta.get('outputs_sha256')
    if not isinstance(outputs, dict):
        raise ValueError(f'Missing stage artifact checksums: {directory}')
    required = ['alignments.bam', 'accepted.bam', 'qc.json']
    if handoff:
        required.append('unmapped.fastq')
    for name in required:
        expected = _digest(outputs.get(name), f'{directory}/{name}')
        path = directory/name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f'Stage artifact checksum mismatch or missing artifact: {path}')
        pins[str(path)] = expected
    qc = _json_object(directory/'qc.json')
    if meta.get('qc') != qc or type(meta.get('input_reads')) is not int or meta['input_reads'] < 0:
        raise ValueError(f'Inconsistent stage QC/provenance: {directory}')
    outcomes = qc.get('primary_outcomes')
    if not isinstance(outcomes, dict) or any(type(n) is not int or n < 0 for n in outcomes.values()):
        raise ValueError(f'Invalid primary-read counts in stage QC: {directory}')
    if meta['input_reads'] != sum(outcomes.values()):
        raise ValueError(f'Inconsistent stage input count: {directory}')
    return meta, qc, original, pins


def verify_handoff_chain(rna, genome):
    """Require RNA output == genome input == current FASTQ byte checksum."""
    rna, genome = Path(rna).resolve(), Path(genome).resolve()
    rm, rq, original_rna, rp = verify_stage(rna, handoff=True)
    gm, gq, _, gp = verify_stage(genome)
    fastq = gm['parameters'].get('fastq')
    inputs = gm.get('inputs_sha256')
    if not isinstance(fastq, str) or not fastq or not isinstance(inputs, dict):
        raise ValueError('Missing genome FASTQ input provenance')
    source = str(Path(fastq).resolve())
    allowed = {str(rna/'unmapped.fastq'), str(original_rna/'unmapped.fastq')}
    if source not in allowed:
        raise ValueError('Genome input path is not the recorded RNA handoff or its verified publication')
    rna_digest = rm['outputs_sha256']['unmapped.fastq']
    genome_digest = _digest(inputs.get(source), 'genome handoff input')
    current_digest = rp[str(rna/'unmapped.fastq')]
    if not rna_digest == genome_digest == current_digest:
        raise ValueError('RNA output / genome input / current handoff checksum mismatch')
    return rq, gq, {**rp, **gp}, {
        'rna_output_sha256': rna_digest, 'genome_input_sha256': genome_digest,
        'current_handoff_sha256': current_digest, 'genome_recorded_fastq': source,
        'stage_provenance_sha256': {'rna': rp[str(rna/'provenance.json')],
                                    'genome': gp[str(genome/'provenance.json')]}}


def accounting(rna, genome, output, log, trim_report=None):
    """Join primary identities on disk; require an exclusive, exhaustive partition."""
    rna, genome, output = Path(rna), Path(genome), Path(output)
    rq, gq, pins, chain = verify_handoff_chain(rna, genome)
    counts = Counter()
    with tempfile.TemporaryDirectory(prefix='account-', dir=output.parent) as tmp:
        with closing(sqlite3.connect(Path(tmp)/'ids.sqlite')) as db:
            db.execute('create table reads(id text primary key, rna text, genome text, accepted integer default 0)')
            for f in bam_records(rna/'alignments.bam', log):
                reason = classify(int(f[1]), int(f[4]), rq['mapq_minimum'], rq['excluded_flags'])
                if reason != 'nonprimary':
                    db.execute('insert into reads(id,rna) values (?,?)', (f[0], reason))
            db.execute('alter table reads add column forwarded integer default 0')
            forwarded = 0
            for name, seq, qual in fastq_records(rna/'unmapped.fastq'):
                row = db.execute('select rna,forwarded from reads where id=?', (name,)).fetchone()
                if row != ('unmapped', 0):
                    raise ValueError('Invalid/duplicate first-stage handoff ID')
                db.execute('update reads set forwarded=1 where id=?', (name,)); forwarded += 1
            for f in bam_records(genome/'alignments.bam', log):
                flag = int(f[1])
                if flag & (256|2048):
                    continue
                row = db.execute('select rna,genome,forwarded from reads where id=?', (f[0],)).fetchone()
                if row != ('unmapped', None, 1):
                    raise ValueError('Genome primary IDs must be exactly first-stage unmapped IDs')
                db.execute('update reads set genome=? where id=?', ('unmapped' if flag & 4 else 'mapped', f[0]))
            missing = db.execute("select count(*) from reads where rna='unmapped' and (genome is null or forwarded=0)").fetchone()[0]
            if missing:
                raise ValueError('Incomplete genome handoff accounting')
            for space, directory in [('rna', rna), ('genome', genome)]:
                n = 0
                for f in bam_records(directory/'accepted.bam', log):
                    row = db.execute('select rna,genome,accepted from reads where id=?', (f[0],)).fetchone()
                    if row is None or row[2] or (space == 'rna' and row[0] != 'accepted') or (space == 'genome' and row[:2] != ('unmapped','mapped')):
                        raise ValueError('Accepted stage identities overlap or violate handoff')
                    db.execute('update reads set accepted=1 where id=?', (f[0],)); n += 1
                if n != (rq if space == 'rna' else gq)['primary_outcomes'].get('accepted', 0):
                    raise ValueError('Accepted BAM/QC mismatch')
                counts[space+'_accepted'] = n
            for reason, n in db.execute('select rna,count(*) from reads group by rna'):
                if reason != 'unmapped':
                    counts['rna_'+reason] = n
            for reason, n in gq['primary_outcomes'].items():
                counts['genome_'+reason] = n
            total = db.execute('select count(*) from reads').fetchone()[0]
            if forwarded != sum(gq['primary_outcomes'].values()) or total != sum(counts.values()):
                raise ValueError('Exclusive stage read accounting failed')
    result = {'status': 'complete', 'trimmed_input_reads': total, 'forwarded_unmapped_reads': forwarded,
              'exclusive_outcomes': dict(counts), 'accepted_read_id_overlap': 0,
              'handoff_sequence_quality_verified': True,
              'verification': 'Complete stage provenance and publication verified; recorded RNA FASTQ output, genome FASTQ input and current handoff SHA-256 agree; consumed BAM/QC hashes verified before and after read-ID accounting',
              'stage_provenance_sha256': chain['stage_provenance_sha256'],
              'handoff_integrity': chain, 'verified_artifacts_sha256': pins}
    if trim_report:
        reads = json.loads(Path(trim_report).read_text())['read_counts']
        if reads['output'] != total or reads['input'] < total:
            raise ValueError('Cutadapt output and RNA-stage read counts differ')
        result.update(raw_input_reads=reads['input'], trimming_removed_reads=reads['input']-total,
                      exclusive_raw_outcomes=dict(counts, trimming_removed=reads['input']-total),
                      trim_report_sha256=sha256(trim_report))
    for path, digest in pins.items():
        if not Path(path).is_file() or sha256(path) != digest:
            raise ValueError(f'Accounting input changed during verification: {path}')
    output.write_text(json.dumps(result, indent=2)+'\n')
    return result
