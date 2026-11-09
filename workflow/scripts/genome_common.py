"""Small shared I/O, safety, QC and provenance helpers for two genome CLIs."""
import argparse
from collections import Counter
from contextlib import closing
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile

from bam_count import command, policy, reference_sizes
from core import sha256

ROOT = Path(__file__).resolve().parents[2]


def base_arguments(p, mapping=False):
    p.add_argument('--output-dir', required=True, help='New external output directory; never overwritten')
    p.add_argument('--threads', type=int, default=1, help='Sequential align/build/sort CPU budget (default 1)')
    p.add_argument('--reference', help='Read-only genome FASTA; required for prepare, optional with existing index')
    p.add_argument('--reference-sha256', help='Optional required FASTA SHA-256')
    p.add_argument('--decoys', help='Exact decoy identifiers, one per line (no name heuristics)')
    p.add_argument('--no-decoys', action='store_true', help='Explicitly declare no decoys')
    if mapping:
        p.add_argument('--fastq', required=True, help='Single-end trimmed 4-line FASTQ, optionally gzip')
        p.add_argument('--index', help='Existing complete read-only index prefix')
        p.add_argument('--mapq', type=int, default=20)
        p.add_argument('--exclude-duplicates', action='store_true')
        p.add_argument('--sort-memory', default='256M', help='SAMtools memory per sort worker (default 256M)')
        p.add_argument('--uncertain-cigar', choices=['reject', 'error'], default='reject')
        p.add_argument('--seed', type=int, default=0, help='Recorded deterministic aligner seed')


def new_output(a):
    if a.threads < 1 or a.threads > 1024:
        raise ValueError('threads must be positive and <=1024')
    out = Path(a.output_dir).expanduser().resolve()
    if out == ROOT or ROOT in out.parents or out.exists():
        raise ValueError('output-dir must be new and outside the source checkout')
    sources = [getattr(a, k, None) for k in ('fastq', 'reference', 'decoys', 'gtf', 'splice_sites', 'exons')]
    if any(out == Path(s).resolve() or out in Path(s).resolve().parents for s in sources if s):
        raise ValueError('Inputs cannot be inside output-dir')
    if a.decoys and a.no_decoys:
        raise ValueError('Choose --decoys or --no-decoys, not both')
    if getattr(a, 'reference_sha256', None) and not a.reference:
        raise ValueError('reference-sha256 requires reference')
    out.mkdir(parents=True)
    return out


def pin(paths):
    return {str(Path(p).resolve()): sha256(p) for p in paths if p}


def check_pins(before):
    if any(sha256(p) != h for p, h in before.items()):
        raise ValueError('Read-only input changed during execution')


def executable(tool):
    path = shutil.which(tool)
    if not path:
        raise ValueError('Required executable missing: '+tool)
    p = subprocess.run([tool, '--version'], text=True, capture_output=True)
    return {'path': str(Path(path).resolve()), 'sha256': sha256(Path(path).resolve()),
            'version_stdout': p.stdout.strip(), 'version_stderr': p.stderr.strip(),
            'version_exit_status': p.returncode}


def index_files(prefix, kind):
    endings = ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2') if kind == 'bowtie2' else tuple('.'+str(i) for i in range(1,9))
    formats = ('bt2', 'bt2l') if kind == 'bowtie2' else ('ht2', 'ht2l')
    families = [[Path(str(prefix)+e+'.'+f).resolve() for e in endings] for f in formats]
    present = [f for f in families if any(p.exists() for p in f)]
    if len(present) != 1 or not all(p.is_file() for p in present[0]):
        raise ValueError('Require one complete unambiguous '+kind+' index family')
    return present[0]


def inspector_tool(prefix, kind):
    if kind=='hisat2':
        native='hisat2-inspect-'+('l' if str(index_files(prefix,kind)[0]).endswith('.ht2l') else 's')
        if shutil.which(native):return native
    return kind+'-inspect'


def dictionary(prefix, kind, out, log, reference=None, reference_hash=None):
    # Prefer the matching native HISAT2 inspector; some vendor Python wrappers
    # still import imp, which is unavailable in Python >=3.12. No vendor edits.
    inspector = inspector_tool(prefix, kind)
    summary = command([inspector, '-s', str(prefix)], log, stdout=subprocess.PIPE, text=True).stdout
    (out/'index.summary.txt').write_text(summary)
    sizes = {}
    for line in summary.splitlines():
        if line.startswith('Sequence-'):
            fields = line.split('\t')
            # Inspectors preserve FASTA descriptions; SAM and FASTA dictionary
            # identifiers use only the first whitespace-delimited token.
            name, length = fields[-2].split()[0], int(fields[-1])
            if name in sizes or length < 1:
                raise ValueError('Invalid index dictionary')
            sizes[name] = length
    if not sizes:
        raise ValueError('Inspector returned no reference dictionary')
    if reference:
        expected = reference_sizes(reference, out/'fasta.sizes.tsv', reference_hash)
        if sizes != expected:
            raise ValueError('FASTA/index identifier or length mismatch')
    (out/'reference.sizes.tsv').write_text(''.join(f'{n}\t{v}\n' for n,v in sizes.items()))
    return sizes


def decoy_names(a, sizes, required=False):
    if required and not a.decoys and not a.no_decoys:
        raise ValueError('Explicit --decoys or --no-decoys required for HISAT2')
    names = Path(a.decoys).read_text().splitlines() if a.decoys else []
    if any(not n or n != n.strip() for n in names) or len(set(names)) != len(names) or not set(names) <= set(sizes):
        raise ValueError('Decoy list must contain unique exact dictionary identifiers')
    return set(names)


def fastq_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as f:
        while True:
            h = f.readline()
            if not h:
                return
            seq, sep, qual = f.readline().rstrip('\r\n'), f.readline().rstrip('\r\n'), f.readline().rstrip('\r\n')
            if not h.startswith('@') or not h[1:].strip() or not sep.startswith('+') or not seq or len(seq)!=len(qual):
                raise ValueError('Malformed 4-line FASTQ')
            if not re.fullmatch('[ACGTNacgtn]+', seq) or any(not 33<=ord(q)<=126 for q in qual):
                raise ValueError('Invalid FASTQ sequence/Phred33 quality')
            yield h[1:].split()[0], seq, qual


def fingerprint(seq, quality):
    return hashlib.sha256((seq+'\n'+quality).encode('ascii')).hexdigest()


def ingest_fastq(source, db):
    db.execute('create table reads (id text primary key, fingerprint text, seen integer default 0, unmapped integer default 0, forwarded integer default 0)')
    n = 0
    try:
        for name, seq, quality in fastq_records(source):
            db.execute('insert into reads (id,fingerprint) values (?,?)', (name, fingerprint(seq, quality))); n+=1
    except sqlite3.IntegrityError as e:
        raise ValueError('FASTQ names must be unique (first whitespace-delimited token)') from e
    db.commit()
    return n


def certain_cigar(cigar, pos, sizes, chrom, seq):
    pieces = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if ''.join(n+op for n,op in pieces)!=cigar or not pieces or any(int(n)<1 for n,op in pieces):
        raise ValueError('Malformed mapped CIGAR')
    ops = [op for n,op in pieces]
    span = sum(int(n) for n,op in pieces if op in 'MDN=X')
    if chrom not in sizes or pos<1 or span<1 or pos+span-1>sizes[chrom]:
        raise ValueError('Alignment outside reference bounds')
    if seq!='*' and sum(int(n) for n,op in pieces if op in 'MIS=X')!=len(seq):
        raise ValueError('CIGAR/query length mismatch')
    return ops[0] in 'M=X' and ops[-1] in 'M=X' and not any(op in 'SHP' for op in ops)


def process_sam(sam, accepted_sam, sizes, decoys, a, db):
    flags = policy(a.mapq, a.exclude_duplicates)
    outcomes, per_reference = Counter(), Counter()
    records = Counter()
    with open(sam) as inp, open(accepted_sam, 'w') as dest:
        for line in inp:
            if line.startswith('@'):
                dest.write(line); continue
            f = line.rstrip('\n').split('\t')
            if len(f)<11:
                raise ValueError('Malformed SAM record')
            name, flag, chrom, pos, mapq = f[0], int(f[1]), f[2], int(f[3]), int(f[4])
            if flag&1:
                raise ValueError('Paired records are outside single-end interface')
            if not 0<=mapq<=255:
                raise ValueError('Invalid SAM MAPQ')
            records['all']+=1
            if flag & (256|2048):
                records['nonprimary']+=1; continue
            row = db.execute('select seen from reads where id=?',(name,)).fetchone()
            if row is None or row[0]:
                raise ValueError('Missing/duplicate primary SAM read identity: '+name)
            db.execute('update reads set seen=1,unmapped=? where id=?',(int(bool(flag&4)),name))
            if flag&4:
                reason='unmapped'
            else:
                per_reference[chrom]+=1
                certain = certain_cigar(f[5],pos,sizes,chrom,f[9])
                if not certain and a.uncertain_cigar=='error':
                    raise ValueError('Coordinate-uncertain CIGAR: '+name+' '+f[5])
                if flag&flags: reason='flags'
                elif chrom in decoys: reason='decoy'
                elif not certain: reason='uncertain_cigar'
                elif mapq<a.mapq: reason='low_mapq'
                else: reason='accepted'
            outcomes[reason]+=1
            if reason=='accepted': dest.write(line)
    if db.execute('select count(*) from reads where seen=0').fetchone()[0]:
        raise ValueError('Aligner omitted input reads')
    return {'primary_outcomes':dict(outcomes),'alignment_records':dict(records),
            'mapped_primary_by_reference':dict(per_reference),
            'decoy_mapped_primary':sum(per_reference[n] for n in decoys),
            'mapq_minimum':a.mapq,'excluded_flags':flags,'decoys':sorted(decoys),
            'cigar_policy':a.uncertain_cigar,'coordinate_policy':'terminal M/= /X required; no S/H/P; internal I/D/N allowed'}


def verify_handoff(path, db):
    count = 0
    for name, seq, quality in fastq_records(path):
        row = db.execute('select fingerprint,unmapped,forwarded from reads where id=?',(name,)).fetchone()
        if row is None or row!=(fingerprint(seq,quality),1,0):
            raise ValueError('Unmapped handoff must be exact original aligner-unmapped reads only')
        db.execute('update reads set forwarded=1 where id=?',(name,)); count+=1
    if count!=db.execute('select count(*) from reads where unmapped=1').fetchone()[0]:
        raise ValueError('Incomplete aligner-unmapped handoff')
    return count


def map_run(a, kind, align, build=None):
    policy(a.mapq,a.exclude_duplicates)
    if not re.fullmatch(r'[1-9][0-9]*[KMG]',a.sort_memory):
        raise ValueError('sort-memory must be a positive SAMtools K/M/G value')
    out = new_output(a)
    paths = [a.fastq,a.reference,a.decoys]
    before = pin(paths)
    meta={'status':'running','tool':kind,'operation':'map','parameters':vars(a),'inputs_sha256':before,
          'python':sys.version,'code_sha256':{p.name:sha256(p) for p in (Path(__file__),ROOT/'workflow/scripts'/f'{kind}_genome.py',ROOT/'workflow/scripts/core.py',ROOT/'workflow/scripts/bam_count.py')},
          'versions':{},'reference_space':'genome; do not combine coordinate keys with collapsed ncRNA matrices'}
    try:
        with (out/'commands.jsonl').open('w') as log, tempfile.TemporaryDirectory(prefix='map-',dir=out) as tmp:
            tmp = Path(tmp)
            # sqlite3's transaction context does not close the connection.
            # Close before removing temp files, including on NFS.
            with closing(sqlite3.connect(tmp/'reads.sqlite')) as db:
                count = ingest_fastq(a.fastq,db); meta['input_reads']=count
                if not a.index:
                    if not build or not a.reference: raise ValueError('Existing index required')
                    a.index=str(out/'index/genome');Path(a.index).parent.mkdir();build(a,log)
                files=index_files(a.index,kind);before.update(pin(files));meta['inputs_sha256']=dict(before)
                sizes=dictionary(a.index,kind,out,log,a.reference,a.reference_sha256)
                decoys=decoy_names(a,sizes,kind=='hisat2')
                native=kind+'-align-'+('l' if str(files[0]).endswith(('bt2l','ht2l')) else 's')
                for tool in [kind,inspector_tool(a.index,kind),'samtools']+([native] if shutil.which(native) else [])+([kind+'-build'] if build and str(out/'index') in a.index else []):
                    meta['versions'][tool]=executable(tool)
                sam=tmp/'all.sam'; accepted=tmp/'accepted.sam'
                align(a,sam,out/'unmapped.fastq',log)
                qc=process_sam(sam,accepted,sizes,decoys,a,db)
                qc['forwarded_unmapped_reads']=verify_handoff(out/'unmapped.fastq',db)
                assert sum(qc['primary_outcomes'].values())==count
                (out/'qc.json').write_text(json.dumps(qc,indent=2)+'\n')
                with (out/'mapping_by_reference.tsv').open('w') as f:
                    f.write('reference\tlength\tdecoy\tmapped_primary\n')
                    for name,length in sizes.items():f.write(f'{name}\t{length}\t{int(name in decoys)}\t{qc["mapped_primary_by_reference"].get(name,0)}\n')
                for source,name in [(sam,'alignments'),(accepted,'accepted')]:
                    # Sequential sort; at most threads total including main SAMtools thread.
                    command(['samtools','sort','-@',str(a.threads-1),'-m',a.sort_memory,'-T',str(tmp/name),'-o',str(out/(name+'.bam')),str(source)],log)
                    command(['samtools','quickcheck','-v',str(out/(name+'.bam'))],log)
                meta['qc']=qc
                db.commit()
        check_pins(before);meta.update(status='complete',inputs_unchanged=True)
        meta['outputs_sha256']={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()}
    except Exception as e:
        meta.update(status='failed',error=str(e));raise
    finally:
        (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    return out


def prepare_run(a, kind, build):
    if not a.reference: raise ValueError('prepare requires --reference')
    out=new_output(a);before=pin([a.reference,a.decoys,getattr(a,'gtf',None),getattr(a,'splice_sites',None),getattr(a,'exons',None)])
    meta={'status':'running','operation':'prepare','tool':kind,'parameters':vars(a),'inputs_sha256':before,'versions':{},
          'python':sys.version,'code_sha256':{p.name:sha256(p) for p in (Path(__file__),ROOT/'workflow/scripts'/f'{kind}_genome.py',ROOT/'workflow/scripts/core.py',ROOT/'workflow/scripts/bam_count.py')}}
    try:
        with (out/'commands.jsonl').open('w') as log:
            sizes=reference_sizes(a.reference,out/'reference.sizes.tsv',a.reference_sha256)
            decoys=decoy_names(a,sizes,kind=='hisat2');a.index=str(out/'index/genome');Path(a.index).parent.mkdir()
            build(a,log)
            dictionary(a.index,kind,out,log,a.reference,a.reference_sha256)
            meta['decoys']=sorted(decoys);meta['index_files_sha256']=pin(index_files(a.index,kind))
            meta['versions']={tool:executable(tool) for tool in [kind+'-build',inspector_tool(a.index,kind)]}
            if getattr(a,'gtf',None):
                meta['annotation_helper_sha256']={str(Path(shutil.which(t)).resolve()):sha256(Path(shutil.which(t)).resolve()) for t in ['hisat2_extract_splice_sites.py','hisat2_extract_exons.py']}
        check_pins(before);meta.update(status='complete',inputs_unchanged=True)
    except Exception as e:
        meta.update(status='failed',error=str(e));raise
    finally:(out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    return out


def cli_entry(parser, runner):
    a=parser.parse_args()
    try:runner(a)
    except (ValueError,OSError,subprocess.SubprocessError) as e:parser.exit(1,f'{parser.prog}: {e}\n')
