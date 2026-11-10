#!/usr/bin/env python3
"""Compact, record-only CapZyme-DAD reports. Standard library; no analysis scans."""
import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
from urllib.parse import quote

VERIFIED = 'Complete and verified (report-time record checks)'
RECORDED = 'Complete according to recorded checks; no fresh data revalidation'
NA = 'Not requested / not applicable'
FAILED = 'Missing, inconsistent or failed'
LIMIT = 16 * 1024 * 1024
METHODS = ('raw', 'rpm', 'deseq2', 'controls', 'supplied')


def workflow_context(settings, config, tracks, outputs, route):
    """Capture requests, not success. No input validation, hashing or tool calls."""
    root = Path(settings['output'])
    return {'schema_version': 1, 'route': route, 'sample_order': list(settings['samples']),
            'settings': dict(settings), 'bigwigs': config.get('bigwigs', {}),
            'requested_outputs': [str(Path(f).relative_to(root)) for f in outputs],
            'track_outputs': {k: [str(Path(f).relative_to(root)) for f in p['files']] for k,p in tracks.items()}}


def records(text):
    """Read JSON objects amid native-tool text/JSON-lines/pretty JSON logs."""
    decoder = json.JSONDecoder(); start = 0
    while start < len(text):
        pos = text.find('{', start)
        if pos < 0: return
        try: obj, end = decoder.raw_decode(text, pos)
        except ValueError: start = pos+1; continue
        start = end
        if isinstance(obj, dict): yield obj


def integer(value):
    if type(value) is not int or value < 0: raise ValueError('Expected nonnegative integer count')
    return value


def digest(value):
    return isinstance(value, str) and bool(re.fullmatch('[0-9a-fA-F]{64}', value))


def encoded(value):
    return html.escape(str(value), quote=True)


class Report:
    def __init__(self, analysis_dir, output_dir=None, context=None):
        self.root = Path(analysis_dir).expanduser().resolve()
        self.out = Path(output_dir).expanduser().resolve() if output_dir else self.root/'report'
        self.checks = []; self.warnings = []; self.links = {}; self.factors = []; self.observed = {}
        self.context = context
        self.names=[]; self.route='unknown'; self.settings={}; self.methods=[]; self.run={}; self.genome={}
        self.rows = []; self.versions = {}; self.references = {}; self.stage_qc = {}

    def check(self, label, ok, detail='', recorded=False):
        self.checks.append({'check': label, 'status': (RECORDED if recorded else VERIFIED) if ok else FAILED,
                            'detail': detail})

    def read(self, path, required=True):
        path = Path(path)
        try:
            if path.stat().st_size > LIMIT: raise ValueError('Record exceeds 16 MiB reporting bound')
            data = path.read_bytes()
            self.observed[str(path)] = hashlib.sha256(data).hexdigest()
            return data.decode('utf-8')
        except (OSError, ValueError) as exc:
            if required: self.check(str(path.relative_to(self.root)) if path.is_relative_to(self.root) else path.name, False, str(exc))
            else: self.warnings.append('Optional record unavailable: '+path.name)
            return None

    def obj(self, path, required=True):
        text = self.read(path, required)
        if text is None: return {}
        try:
            data = json.loads(text)
            if not isinstance(data, dict): raise ValueError('Expected JSON object')
            return data
        except ValueError as exc:
            self.check(Path(path).name, False, str(exc)); return {}

    def output(self, relative, explanation):
        p = self.root/relative
        # Only inspect paths/headers; never read BAM/FASTQ/BigWig bodies.
        ok = p.is_file()
        self.check('Output '+relative, ok, 'Present; contents/hashes not freshly validated' if ok else 'Required output absent')
        self.links[relative] = explanation
        return ok

    def log_result(self, path, step=None):
        text = self.read(path)
        if text is None: return {}
        found = [o for o in records(text) if 'exit_status' in o and (step is None or o.get('step') == step)]
        if not found or found[-1].get('exit_status') != 0:
            self.check(Path(path).name, False, 'Missing/failed completion record'); return {}
        self.check(Path(path).name, True, 'Saved step exit status is zero', recorded=True)
        return found[-1].get('summary', found[-1].get('result', {}))

    def matrix(self, relative):
        if not self.output(relative, 'Integer one-based 5′-end counts; statistical input, not normalized tracks'): return
        try:
            with (self.root/relative).open() as f: header = f.readline(65537).rstrip('\r\n').split('\t')
            self.check(relative+' sample columns', header == ['chromosome','strand','start','end',*self.names],
                       'Fresh header check only; feature rows/count cells not rescanned')
        except (OSError, UnicodeError) as exc: self.check(relative, False, str(exc))

    def load(self):
        if self.context is None:
            saved = self.root/'report/analysis_context.json'
            self.context = self.obj(saved) if saved.exists() else {}
        c = self.context
        if not isinstance(c,dict):
            self.check('Request plan',False,'Expected context object'); c={}; self.context=c
        run_path=self.root/'provenance/run.json'
        self.run = self.obj(run_path) if run_path.exists() or c.get('route')!='genome' and not (self.root/'provenance/genome.json').exists() else {}
        self.genome = self.obj(self.root/'provenance/genome.json', required=c.get('route') in ('genome','sequential')) if (self.root/'provenance/genome.json').exists() or c.get('route') in ('genome','sequential') else {}
        s = self.run.get('settings') or self.genome.get('settings') or c.get('settings') or {}
        self.settings = s
        if c:
            self.check('Request plan schema',type(c.get('schema_version')) is int and c.get('schema_version')==1,'Saved request plan is not a workflow completion receipt')
        if c.get('settings') and (self.run or self.genome):
            fields=('pre','mapq','flags','strand','minimum','genome','reference_mode','fixed','fasta')
            self.check('Requested versus recorded settings', all(c['settings'].get(k)==s.get(k) for k in fields), 'Analysis settings compared; track requests may be added independently')
        inferred = ('sequential' if s.get('genome') and 'pre' in s else 'genome' if s.get('genome') else 'rna' if 'pre' in s else 'bam' if s.get('fasta') else None)
        self.route = c.get('route') or inferred or 'unknown'
        self.names = c.get('sample_order') or self.run.get('sample_order') or list(s.get('samples', {}))
        self.check('Ordered sample IDs', isinstance(self.names,list) and bool(self.names) and all(isinstance(n,str) and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',n) for n in self.names) and len(set(self.names))==len(self.names), 'IDs retained as strings; no treatment inferred')
        if not isinstance(self.names,list) or any(not isinstance(n,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',n) for n in self.names): self.names=[]
        if inferred: self.check('Route identity', self.route == inferred, 'Requested versus recorded route')
        for label, obj in [('RNA/run provenance',self.run),('Genome provenance',self.genome)]:
            if not obj: continue
            order = obj.get('sample_order') or list(obj.get('settings',{}).get('samples',{}))
            self.check(label+' sample order', order == self.names, 'Recorded order versus report request')
            source = obj.get('git') or obj.get('source_snapshot') or obj.get('installed_source') or {}
            self.versions[label] = {'source': {k: source[k] for k in ('commit','branch','dirty','version','distribution','kind') if k in source},
                'snapshot_sha256':obj.get('source_snapshot_sha256',obj.get('execution_snapshot_sha256')),
                'software':obj.get('versions', {k:obj[k] for k in ('python','snakemake') if k in obj})}
        if not c:
            self.warnings.append('No saved request plan: route inferred from provenance; enabled track methods discovered from saved publications. Undiscovered historical requests cannot be certified.')
        if self.route!='sequential':
            self.checks.append({'check':'RNA/genome overlap and handoff','status':NA,'detail':'No sequential handoff requested for '+self.route})
        self.check('Analysis route available', self.route in ('rna','sequential','genome','bam'), self.route)
        for label,receipt in [('run',self.run)] if self.route!='genome' else [('genome',self.genome)]:
            self.check(label+' provenance structure',isinstance(receipt.get('settings'),dict) and list(receipt.get('settings',{}).get('samples',{}))==self.names and bool(self.names),'Required saved cohort settings; an existing file alone is insufficient')
        self.metadata = {}
        manifest = s.get('manifest')
        if manifest and Path(manifest).is_file():
            text = self.read(manifest, required=False)
            if text:
                rows = list(csv.DictReader(text.splitlines(), delimiter='\t'))
                self.check('Manifest sample order', [r.get('sample_id') for r in rows] == self.names, 'Existing manifest strings; no input data opened')
                self.metadata = {r['sample_id']:{k:r[k] for k in ('treatment','replicate') if r.get(k)} for r in rows if 'sample_id' in r}
        for space in ('rna','genome'):
            ref = s.get('genome',{}) if space=='genome' else s.get('fixed') or {'fasta':s.get('fasta'), 'sources':s.get('reference_sources')}
            if ref: self.references[space] = {k:v for k,v in ref.items() if k in ('fasta','reference','reference_sha256','fasta_sha256','index_prefix','index_format','sources','policy_files','decoys')}
        self.references['recorded_input_sha256'] = {k:v for k,v in self.run.get('inputs',{}).items() if k in s.get('reference_sources',[])}

    def stage(self, space, name):
        d = self.root/space/name; q = self.obj(d/'qc.json'); m = self.obj(d/'provenance.json')
        self.output(f'{space}/{name}/accepted.bam', 'Accepted single-end alignments; use matching reference')
        self.links[f'{space}/{name}/qc.json']='Detailed stage accounting and filters'
        self.links[f'{space}/{name}/commands.jsonl']='Native commands, versions, exit statuses and logs'
        try:
            outcomes = q['primary_outcomes']; n = integer(m['input_reads'])
            if not isinstance(outcomes,dict): raise ValueError('Missing primary outcome map')
            for k,v in outcomes.items(): integer(v)
            if sum(outcomes.values()) != n: raise ValueError('Stage outcomes do not sum to stage input')
            self.check(name+' '+space+' accounting', True, 'Fresh sum of recorded exclusive primary outcomes')
            self.check(name+' '+space+' completion', m.get('status')=='complete' and m.get('inputs_unchanged') is True and m.get('qc')==q,
                       'Saved completion/input verification and QC identity', recorded=True)
            outputs = m.get('outputs_sha256',{})
            self.check(name+' '+space+' QC checksum', outputs.get('qc.json') == self.observed.get(str(d/'qc.json')), 'Fresh small QC-file hash versus stage receipt')
            for filename,pin in outputs.items():
                if not isinstance(filename,str) or Path(filename).is_absolute() or '..' in Path(filename).parts:
                    self.check(name+' '+space+' output pin path',False,'Expected relative artifact name'); continue
                self.output(f'{space}/{name}/{filename}','Stage artifact recorded by completed mapper; body not reread')
                self.check(name+' '+space+' pin '+filename,digest(pin),'Saved SHA-256 syntax only',recorded=True)
            self.check(name+' '+space+' BAM pins recorded', all(digest(outputs.get(k)) for k in ('accepted.bam','alignments.bam')), 'Saved hashes only; BAMs not read', recorded=True)
            pub = self.obj(d/'publication.json')
            self.check(name+' '+space+' publication binding', pub.get('standalone_provenance_sha256')==self.observed.get(str(d/'provenance.json')) and pub.get('standalone_output_dir')==m.get('parameters',{}).get('output_dir') and pub.get('published_output_dir')==str(d), 'Fresh provenance hash and recorded publication locations')
            versions = m.get('versions',{})
            self.versions[space] = {k: (v.get('version_stdout','').splitlines() or ['unavailable'])[0] if isinstance(v,dict) else v for k,v in versions.items()}
            pars=m.get('parameters',{}); ref=pars.get('reference'); self.references[space]={'reference':ref,'recorded_sha256':m.get('inputs_sha256',{}).get(ref),'index':pars.get('index'),**self.references.get(space,{})}
            nonprimary=integer(q.get('alignment_records',{}).get('nonprimary',0))
            all_records=q.get('alignment_records',{}).get('all')
            if all_records is not None: self.check(name+' '+space+' alignment records', integer(all_records)==n+nonprimary, 'Nonprimary alignment records are separate from primary input-read totals')
            self.stage_qc[space,name]=(q,m)
            return n,outcomes,nonprimary
        except (KeyError,ValueError,TypeError,AttributeError) as exc:
            self.check(name+' '+space+' accounting',False,str(exc));return None,{},None

    def libraries(self):
        for name in self.names:
            row={'sample_id':name,**self.metadata.get(name,{})}
            for k in ('raw_reads','trimming_removed_reads','trimmed_reads','rna_input_reads','rna_mapped_reads','rna_accepted_reads','rna_low_mapq_reads','rna_flag_rejected_reads','rna_unmapped_reads','rna_nonprimary_records','genome_input_reads','genome_mapped_reads','genome_accepted_reads','genome_low_mapq_reads','genome_flag_rejected_reads','genome_decoy_rejected_reads','genome_cigar_rejected_reads','genome_nonprimary_records','final_unmapped_reads','rna_counted_ends','genome_counted_ends'): row[k]=None
            if self.route in ('rna','sequential'):
                trim=self.obj(self.root/'qc'/f'{name}.cutadapt.json')
                try:
                    rc=trim['read_counts'];raw=integer(rc['input']);out=integer(rc['output'])
                    if raw<out: raise ValueError('Cutadapt output exceeds input')
                    row.update(raw_reads=raw,trimmed_reads=out,trimming_removed_reads=raw-out)
                    self.check(name+' trimming',True,'Fresh input/output difference; counts from saved Cutadapt JSON')
                except (KeyError,ValueError,TypeError,AttributeError) as exc: self.check(name+' trimming',False,str(exc))
            for space in ('rna','genome'):
                if space=='rna' and self.route!='sequential' or space=='genome' and self.route not in ('sequential','genome'): continue
                n,o,np=self.stage(space,name);row[space+'_input_reads']=n;row[space+'_nonprimary_records']=np
                if n is None: continue
                row.update({space+'_mapped_reads':n-o.get('unmapped',0),space+'_accepted_reads':o.get('accepted',0),space+'_low_mapq_reads':o.get('low_mapq',0),space+'_flag_rejected_reads':o.get('flags',0)})
                if space=='rna': row['rna_unmapped_reads']=o.get('unmapped',0)
                else: row.update(final_unmapped_reads=o.get('unmapped',0),genome_decoy_rejected_reads=o.get('decoy',0),genome_cigar_rejected_reads=o.get('uncertain_cigar',0))
            if self.route=='rna':
                self.output(f'bam/{name}.accepted.bam','Accepted RNA alignments; matching RNA reference')
                text=self.read(self.root/'logs'/f'{name}.bowtie2.log')
                if text:
                    match=re.search(r'(\d+) reads; of these:',text);zero=re.search(r'(\d+) \([^\n]*\) aligned 0 times',text)
                    if match and zero:
                        n=int(match[1]);un=int(zero[1]);row.update(rna_input_reads=n,rna_mapped_reads=n-un,rna_unmapped_reads=un,final_unmapped_reads=un)
                        self.check(name+' RNA aligner input', n==row['trimmed_reads'], 'Native Bowtie2 summary versus Cutadapt output; post-filter rejects unavailable')
                self.warnings.append(name+': RNA-only streaming route lacks per-category filtering receipts; unavailable fields remain blank.')
            if self.route in ('rna','sequential') and row['rna_input_reads'] is not None:
                self.check(name+' trim → RNA',row['trimmed_reads']==row['rna_input_reads'],'Recorded stage counts agree')
            if self.route=='bam':
                self.output(f'bam/{name}.accepted.bam','Output-local prepared BAM; original source BAM unchanged')
                prep=self.log_result(self.root/'logs'/f'{name}.prepare.jsonl','prepare')
                self.check(name+' source BAM preservation',prep.get('source_unchanged_stat') is True,'Saved stat check, not fresh source revalidation',recorded=True)
            if self.route=='sequential': self.handoff(name,row)
            self.rows.append(row)
        if self.route=='genome': self.warnings.append('Genome-only input is already trimmed: raw counts and trimming losses are not applicable/unavailable; percentages use genome-stage input.')
        if self.route=='bam': self.warnings.append('Existing-BAM route: FASTQ trimming/mapping accounting is not applicable. Count totals refer to eligible aligned ends on configured strands.')

    def handoff(self,name,row):
        a=self.obj(self.root/'qc'/f'{name}.sequential.json')
        try:
            outcomes=a['exclusive_outcomes'];total=integer(a['trimmed_input_reads'])
            self.check(name+' exclusive accounting', all(type(v) is int and v>=0 for v in outcomes.values()) and sum(outcomes.values())==total==row['trimmed_reads'], 'Fresh reconciliation of saved primary-outcome counts')
            want={k:int(v or 0) for k,v in [('rna_accepted',row['rna_accepted_reads']),('rna_low_mapq',row['rna_low_mapq_reads']),('rna_flags',row['rna_flag_rejected_reads']),('genome_accepted',row['genome_accepted_reads']),('genome_low_mapq',row['genome_low_mapq_reads']),('genome_flags',row['genome_flag_rejected_reads']),('genome_decoy',row['genome_decoy_rejected_reads']),('genome_uncertain_cigar',row['genome_cigar_rejected_reads']),('genome_unmapped',row['final_unmapped_reads'])]}
            self.check(name+' outcome identity', {k:v for k,v in outcomes.items() if v}=={k:v for k,v in want.items() if v}, 'Exclusive accounting versus each stage QC')
            self.check(name+' raw accounting', a.get('raw_input_reads')==row['raw_reads'] and a.get('trimming_removed_reads')==row['trimming_removed_reads'] and sum(a.get('exclusive_raw_outcomes',{}).values())==row['raw_reads'], 'Saved raw partition versus Cutadapt; raw denominator includes trimming losses')
            self.check(name+' forwarded counts', a.get('forwarded_unmapped_reads')==row['rna_unmapped_reads']==row['genome_input_reads'], 'RNA aligner-unmapped only → genome-stage input')
            chain=a['handoff_integrity'];values=[chain.get(k) for k in ('rna_output_sha256','genome_input_sha256','current_handoff_sha256')]
            self.check(name+' saved handoff checksums',all(digest(v) for v in values) and len(set(values))==1,'Compare saved RNA output/genome input/analysis-time current handoff pins; no FASTQ hash read now')
            rm=self.stage_qc.get(('rna',name),({},{}))[1];gm=self.stage_qc.get(('genome',name),({},{}))[1]
            source=gm.get('parameters',{}).get('fastq')
            self.check(name+' handoff receipt binding',values[0]==rm.get('outputs_sha256',{}).get('unmapped.fastq') and values[1]==gm.get('inputs_sha256',{}).get(source) and chain.get('stage_provenance_sha256')==a.get('stage_provenance_sha256') and all(a.get('stage_provenance_sha256',{}).get(space)==self.observed.get(str(self.root/space/name/'provenance.json')) for space in ('rna','genome')), 'Fresh small-receipt hashes match saved accounting chain')
            self.check(name+' accepted read-ID overlap / handoff',a.get('status')=='complete' and a.get('accepted_read_id_overlap')==0 and type(a.get('accepted_read_id_overlap')) is int and a.get('handoff_sequence_quality_verified') is True, 'Recorded overlap=0 and sequence/quality verification; read identities NOT rechecked now',recorded=True)
            self.check(name+' trimming receipt binding',a.get('trim_report_sha256')==self.observed.get(str(self.root/'qc'/f'{name}.cutadapt.json')),'Fresh small Cutadapt JSON hash')
        except (KeyError,ValueError,TypeError,AttributeError) as exc: self.check(name+' sequential integrity',False,str(exc))

    def counts(self):
        for space in ('rna','genome'):
            if space=='rna' and self.route not in ('rna','sequential','bam') or space=='genome' and self.route not in ('genome','sequential'):continue
            relative='counts/'+('genome.5ends.raw.tsv' if space=='genome' else '5ends.raw.tsv');self.matrix(relative)
            path=self.root/'logs'/('count.jsonl' if self.route=='bam' else 'count_5ends.json' if space=='rna' else 'genome.count.jsonl')
            # Genome counter currently records commands/statuses but no cohort summary.
            text=self.read(path)
            if text:
                exits=[int(v) for v in re.findall(r'exit_status=(\d+)',text)]
                json_exits=[o['exit_status'] for o in records(text) if 'exit_status' in o]
                self.check(space+' counter completion records',bool(exits or json_exits) and all(v==0 for v in exits+json_exits),'Saved native subprocess/step exit statuses; not a fresh count scan',recorded=True)
                summaries=[o.get('summary',o.get('result',{})) for o in records(text) if 'accepted_ends_per_sample' in o.get('summary',o.get('result',{}))]
                if summaries:
                    x=summaries[-1];vals=x['accepted_ends_per_sample']
                    ok=x.get('samples')==self.names and len(vals)==len(self.names) and all(type(v) is int and v>=0 for v in vals)
                    self.check(space+' count summary',ok,'Recorded counts and sample order; not a fresh matrix reconstruction')
                    if ok:
                        for row,n in zip(self.rows,vals):
                            row[space+'_counted_ends']=n
                            if self.settings.get('strand')=='both' and space=='rna':
                                if row['rna_accepted_reads'] is None:row['rna_accepted_reads']=n
                                else:self.check(row['sample_id']+' RNA counted totals',n==row['rna_accepted_reads'],'Both-strand saved ends versus saved accepted count')
                elif space=='genome':self.warnings.append('Genome count log has no cohort-total summary; matrix header inspected, counts not reconstructed.')
                else:self.check(space+' count completion',False,'Missing count summary')
            if self.route=='rna':self.log_result(self.root/'logs/provenance.json','provenance')
            if self.route=='bam':self.log_result(self.root/'logs/provenance.jsonl','provenance')
        minimum=self.settings.get('minimum',0)
        if minimum>0 or self.route=='bam':self.matrix('counts/5ends.min_total_'+str(minimum)+'.tsv')
        if self.settings.get('comparison'):
            x=self.obj(self.root/'comparison/exact.json');self.check('Historical comparison',x.get('equal') is True,'Saved keyed comparison; not repeated',recorded=True)

    def normalization(self):
        blocks=self.context.get('bigwigs',{})
        requested={space+'/'+method:options for space,b in blocks.items() for method,options in b.get('methods',{}).items()}
        if not self.context:
            requested={space+'/'+m:{} for space in ('rna','genome') for m in METHODS if (self.root/'bigwigs'/space/m).exists()}
        self.methods=[]; identities={}; sample_hashes={}
        for space in ('rna','genome'):
            for method in METHODS:
                key=space+'/'+method
                if key not in requested:
                    self.methods.append((key,'not enabled' if self.context else 'not enabled in available records (request plan unavailable)'));continue
                self.methods.append((key,'enabled'+(' — external declared depth pool' if requested[key].get('depths') else '')))
                d=self.root/'bigwigs'/space/method;m=self.obj(d/'provenance.json');opts=requested[key]
                try:
                    cohort=m['cohort'];emitted=m['emitted_samples'];pars=m['parameters'];block=blocks.get(space,{})
                    identity=(m.get('matrix_sha256'),m.get('reference_dictionary_sha256'))
                    self.check(key+' normalization source identity',all(digest(x) for x in identity) and (space not in identities or identity==identities[space]),'Saved matrix/reference identities agree across methods within this space; no fresh body hashes')
                    identities.setdefault(space,identity)
                    for sample,pin in m.get('sample_count_hashes',{}).items():
                        self.check(key+' '+sample+' saved count identity',digest(pin) and ((space,sample) not in sample_hashes or sample_hashes[space,sample]==pin),'Cross-method saved count fingerprints only')
                        sample_hashes.setdefault((space,sample),pin)
                    if method in ('deseq2','controls'):
                        self.check(key+' estimator',m.get('estimator') in ('ratio','poscounts') and (not self.context or m.get('estimator')==opts.get('deseq2_type')),'Explicit saved DESeq2 estimator; no R invoked or estimator switched')
                    expected_matrix=self.root/'counts'/('genome.5ends.raw.tsv' if space=='genome' else '5ends.raw.tsv')
                    if pars.get('matrix'):
                        for receipt in (self.run,self.genome):
                            expected_pin=receipt.get('inputs_sha256',{}).get(str(expected_matrix))
                            if expected_pin is not None:self.check(key+' recorded matrix pin',expected_pin==m.get('matrix_sha256'),'Compare two saved records; matrix body not hashed now')
                        self.check(key+' matrix space',Path(pars['matrix']).resolve()==expected_matrix,'Recorded matrix path belongs to this reference space; body hash not rechecked')
                    self.check(key+' completion',m.get('status')=='complete' and m.get('normalization')==method,'Saved normalization completion',recorded=True)
                    ok=isinstance(cohort,list) and bool(cohort) and len(set(cohort))==len(cohort) and all(n in self.names for n in cohort) and isinstance(emitted,list) and bool(emitted) and len(set(emitted))==len(emitted) and all(n in cohort for n in emitted)
                    if self.context:ok=ok and cohort==block.get('cohort',self.names) and emitted==block.get('emit_samples',cohort) and pars.get('strands')==block.get('strands',self.settings.get('genome',{}).get('strand','both') if space=='genome' else self.settings.get('strand','both'))
                    self.check(key+' cohort/selection',ok,'One declared within-space cohort; no RNA/genome pool merging')
                    text=self.read(d/'normalization_factors.tsv')
                    rows=list(csv.DictReader((text or '').splitlines(),delimiter='\t'))
                    self.check(key+' factor identities',[r.get('sample_id') for r in rows]==cohort,'Exact ordered factor cohort, including leading zeros')
                    self.check(key+' factor checksum',m.get('outputs_sha256',{}).get('normalization_factors.tsv')==self.observed.get(str(d/'normalization_factors.tsv')),'Fresh small TSV hash against normalization provenance')
                    if (d/'workflow.json').exists():
                        w=self.obj(d/'workflow.json');self.check(key+' workflow receipt',w.get('status')=='complete' and w.get('exit_status')==0,'Saved standalone BigWig subprocess completion',recorded=True)
                    for r in rows:
                        v=float(r['multiplier']);sf=r.get('size_factor');depth=r.get('accepted_records')
                        valid=math.isfinite(v) and v>0 and r.get('method')==method and r.get('cohort_id')==m.get('cohort_id') and r.get('units')==m.get('units')
                        if method=='raw':valid=valid and v==1
                        if method=='rpm':
                            valid=valid and depth.isdigit() and int(depth)>0 and math.isclose(v,1e6/int(depth),rel_tol=1e-12)
                            saved_row=next((x for x in self.rows if x['sample_id']==r['sample_id']),{})
                            accepted=saved_row.get(space+'_accepted_reads')
                            if accepted is not None and not opts.get('depths'):
                                self.check(key+' '+r['sample_id']+' accepted pool',int(depth)==accepted, 'Saved RPM depth versus saved accepted alignments in this stage, both strands')
                        if method in ('deseq2','controls') or sf:
                            factor=float(sf);valid=valid and math.isfinite(factor) and factor>0 and math.isclose(v,1/factor,rel_tol=1e-12)
                        self.check(key+' '+r['sample_id']+' factor',valid,'Fresh arithmetic/identity checks on saved factors; no normalization recomputed')
                        self.factors.append({'space':space,**r,'denominator':m.get('denominator_definition','not RPM'),'estimator':m.get('estimator','not applicable')})
                    for name in emitted:
                        for strand in ('plus','minus') if pars['strands']=='both' else (pars['strands'],):
                            rel=f'bigwigs/{key}/tracks/{name}.{strand}.bw';self.output(rel,'Single-base reference-relative 5′ signal; '+str(m.get('units')))
                            self.check(rel+' recorded checksum',digest(m.get('outputs_sha256',{}).get('tracks/'+name+'.'+strand+'.bw')),'Saved pin only; BigWig body not decoded now',recorded=True)
                    self.warnings.extend(str(w) for w in m.get('warnings',[]))
                    self.links[f'bigwigs/{key}/normalization_factors.tsv']='Normalization cohort, factors, denominators and units'
                    self.links[f'bigwigs/{key}/provenance.json']='Detailed normalization/source/reference records'
                except (KeyError,ValueError,TypeError,AttributeError) as exc:self.check(key+' normalization',False,str(exc))

    def generate(self):
        for phase in (self.load,self.libraries,self.counts,self.normalization):
            try: phase()
            except (KeyError,ValueError,TypeError,AttributeError,OSError) as exc:
                self.check(phase.__name__+' records',False,str(exc))
        for receipt in (self.run,self.genome):
            for path,pin in receipt.get('inputs_sha256',{}).items():
                if path in self.observed:
                    self.check('Receipt binding '+str(Path(path).name),pin==self.observed[path], 'Fresh bounded record hash; no sequencing/artifact body read')
        for rel in self.context.get('requested_outputs',[]):
            if not Path(rel).is_absolute() and '..' not in Path(rel).parts:self.output(rel,'Requested analysis artifact or completion record')
            else:self.check('Requested output path',False,'Unsafe/nonrelative request path')
        for name in self.names:
            for key in ('rna','genome'):
                row=next(r for r in self.rows if r['sample_id']==name)
                for numerator,denominator in [(key+'_accepted_reads','raw_reads'),(key+'_accepted_reads',key+'_input_reads'),(key+'_mapped_reads',key+'_input_reads')]:
                    num,den=row.get(numerator),row.get(denominator)
                    row[numerator+'_pct_of_'+denominator]=None if num is None or not den else round(100*num/den,6)
        self.out.mkdir(parents=True,exist_ok=True)
        cols=['sample_id']+[k for k in ('treatment','replicate') if any(k in r for r in self.rows)]
        cols+=list(dict.fromkeys(k for r in self.rows for k in r if k not in cols))
        with (self.out/'sample_summary.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:('' if v is None else v) for k,v in r.items()} for r in self.rows)
        failed=any(c['status']==FAILED for c in self.checks)
        status={'schema_version':1,'status':'failed' if failed else 'complete_from_records', 'analysis_route':self.route,'sample_order':self.names,
                'verification_scope':'Bounded saved-record checks, requested-file presence and matrix headers only. No fresh large-file hashes/read IDs/count reconstruction/BigWig decoding.',
                'workflow_exit_status':'not assessed; a final workflow task cannot certify subsequent coordinator exit',
                'checks':self.checks,'warnings':list(dict.fromkeys(self.warnings)),'records_sha256':self.observed}
        (self.out/'report_status.json').write_text(json.dumps(status,indent=2)+'\n')
        if self.context:(self.out/'analysis_context.json').write_text(json.dumps(self.context,indent=2)+'\n')
        (self.out/'report.html').write_text(self.render(status),encoding='utf-8')
        return 1 if failed else 0

    def link(self,p,label=None):
        relative=os.path.relpath(p,self.out);return '<a href="'+encoded(quote(relative,safe='/'))+'">'+encoded(label or relative)+'</a>'

    @staticmethod
    def table(headers,rows):
        return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+encoded(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+encoded('—' if v is None else v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'

    def render(self,status):
        title='Missing, inconsistent or failed required checks' if status['status']=='failed' else RECORDED
        overview={'Ordered samples':', '.join(self.names),'Route':self.route,
                  'Recorded source/software':self.versions,'References (recorded identities; not rehashed)':self.references,
                  'Important settings':{k:self.settings[k] for k in ('pre','reference_mode','norc','preset','seed','mapq','flags','minimum','strand','exclude_duplicates') if k in self.settings},
                  'Genome settings':{k:v for k,v in (self.settings.get('genome') or {}).items() if k in ('aligner','seed','mapq','strand','exclude_duplicates','decoys','no_decoys')}}
        sections='<h1>CapZyme-DAD run report</h1><p class="status">'+encoded(title)+'</p>'
        sections+='<p>'+encoded(status['verification_scope'])+' Workflow exit: '+encoded(status['workflow_exit_status'])+'. Timing: not reported; no reliable overall start/end record is required or inferred from mtimes.</p>'
        if status['warnings']:sections+='<aside><strong>Warnings / incomplete information</strong><ul>'+''.join('<li>'+encoded(w)+'</li>' for w in status['warnings'])+'</ul></aside>'
        errors=[x for x in self.checks if x['status']==FAILED]
        if errors:sections+='<aside class="error"><strong>Required checks failed</strong><ul>'+''.join('<li>'+encoded(x['check']+': '+x['detail'])+'</li>' for x in errors)+'</ul></aside>'
        sections+='<h2>Run overview</h2>'+self.table(['Ordered samples','Mapping route'],[[', '.join(self.names),self.route]])
        sources=[]
        for v in self.versions.values():
            if isinstance(v,dict) and isinstance(v.get('source'),dict):
                source=v['source']; label=' '.join(str(source[k]) for k in ('distribution','version','commit') if source.get(k))
                if label and label not in sources:sources.append(label)
        sections+='<p>Analysis source: '+encoded('; '.join(sources) or 'unavailable in saved records')+'.</p>'
        for key,value in overview.items():
            if key in ('Ordered samples','Route'):continue
            sections+='<details><summary>'+encoded(key)+'</summary><pre>'+encoded(json.dumps(value,indent=2,ensure_ascii=False))+'</pre></details>'
        cols=['sample_id','treatment','replicate','raw_reads','trimming_removed_reads','rna_input_reads','rna_mapped_reads','rna_accepted_reads','rna_low_mapq_reads','rna_flag_rejected_reads','rna_unmapped_reads','genome_input_reads','genome_mapped_reads','genome_accepted_reads','genome_low_mapq_reads','genome_flag_rejected_reads','genome_decoy_rejected_reads','genome_cigar_rejected_reads','final_unmapped_reads','rna_accepted_reads_pct_of_raw_reads','rna_accepted_reads_pct_of_rna_input_reads','genome_accepted_reads_pct_of_raw_reads','genome_accepted_reads_pct_of_genome_input_reads']
        cols=[k for k in cols if k=='sample_id' or any(r.get(k) is not None for r in self.rows)]
        sections+='<h2>Per-library accounting</h2><p>Stage applicability: trimming '+('enabled' if self.route in ('rna','sequential') else 'not applicable')+'; RNA mapping '+('enabled' if self.route in ('rna','sequential') else 'not applicable')+'; genome mapping '+('enabled' if self.route in ('genome','sequential') else 'not applicable')+'.</p><p>RNA/genome mapped rates = mapped primary reads / that stage’s input. Accepted stage rate = accepted reads / that stage’s input. Original-input percentage = accepted reads / raw reads (includes trimming losses). All percentage column names in TSV state their denominator; undefined zero-denominator rates and unavailable/not-applicable fields are blank. Nonprimary records are separate from input-read totals. No treatment/replicate metadata is inferred.</p>'+self.link(self.out/'sample_summary.tsv','Download all counts, explicit percentages and metadata (TSV)')+self.table([k.replace('_reads_pct_of_raw_reads',' % raw').replace('_reads_pct_of_rna_input_reads',' % RNA input').replace('_reads_pct_of_genome_input_reads',' % genome input').replace('_reads','').replace('_',' ').replace('rna','RNA').replace('genome','Genome') for k in cols],[[r.get(k) for k in cols] for r in self.rows])
        sections+='<h2>Integrity and completion checks</h2><p>“Complete and verified” refers only to the specific report-time record check, never fresh validation of sequencing data. Saved identity/large-file checks are explicitly labelled recorded. '+encoded(NA)+'.</p>'+'<p>'+str(sum(c['status'] in (VERIFIED,RECORDED) for c in self.checks))+' checks passed; '+str(len(errors))+' required checks failed.</p><details><summary>Show record checks and their verification scope</summary>'+self.table(['Check','Status','Evidence'],[[c['check'],c['status'],c['detail']] for c in self.checks])+'</details>'
        sections+='<h2>Normalization</h2>'+self.table(['Reference space / method','Request state'],self.methods)
        if self.factors:
            keys=['space','sample_id','method','accepted_records','size_factor','multiplier','units','denominator','estimator','cohort_id'];sections+=self.table(keys,[[r.get(k,'') for k in keys] for r in self.factors])
        sections+='<p>RNA and genome are separate normalization pools. One genomic sample factor applies across its strands. DESeq2 values are counts / size factor, not RPM. Statistical analysis uses original integer counts.</p><h2>Output guide</h2>'
        for label,match in [('Count matrices',lambda p:p.startswith('counts/')),('Accepted BAMs',lambda p:p.endswith('accepted.bam')),('BigWig tracks',lambda p:p.endswith('.bw')),('Detailed records',lambda p:not p.startswith('counts/') and not p.endswith(('accepted.bam','.bw')))]:
            links=[(p,why) for p,why in sorted(self.links.items()) if match(p) and (self.root/p).is_file()]
            sections+='<details'+(' open' if label=='Count matrices' else '')+'><summary>'+encoded(label)+' ('+str(len(links))+')</summary><ul>'
            sections+=''.join('<li>'+self.link(self.root/p,p)+' — '+encoded(why)+'</li>' for p,why in links)+'</ul></details>'
        sections+='<ul>'
        for folder in ('qc','logs','provenance'):
            if (self.root/folder).is_dir():sections+='<li>'+self.link(self.root/folder,folder+'/')+' — detailed native logs, QC and saved provenance</li>'
        sections+='</ul><p>MultiQC is optional and external: run it separately on recognised native logs, keeping RNA and genome stages distinguishable. It is not run or required by this report.</p>'
        return '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CapZyme-DAD run report</title><style>body{font:15px/1.5 system-ui,sans-serif;color:#202b35;max-width:1200px;margin:2rem auto;padding:0 1rem}h2{border-bottom:1px solid #ccd5dc;margin-top:2rem}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f6f8;padding:.7rem}table{border-collapse:collapse;font-size:13px;width:100%}td,th{border:1px solid #ccd5dc;padding:.4rem;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#edf2f6;min-width:6rem}.scroll{overflow-x:auto}aside,.status{background:#fff4d6;border-left:4px solid #ad7d00;padding:1rem}.error{background:#ffe9e6;border-color:#a52820}a{color:#165ca3}</style></head><body>'+sections+'</body></html>\n'


def generate_report(analysis_dir, output_dir=None, context=None):
    return Report(analysis_dir,output_dir,context).generate()


def main():
    p=argparse.ArgumentParser(prog='capzyme report',description=__doc__)
    p.add_argument('--analysis-dir',required=True,help='Existing workflow output directory; never scans sequencing files')
    p.add_argument('--output-dir',help='Report directory (default ANALYSIS_DIR/report); records may contain private paths')
    a=p.parse_args()
    try:code=generate_report(a.analysis_dir,a.output_dir)
    except (OSError,ValueError,TypeError) as exc:p.exit(1,'report: '+str(exc)+'\n')
    print('Report written; '+('required checks failed' if code else 'saved-record checks passed')+'; no analysis rerun')
    raise SystemExit(code)

if __name__=='__main__':main()
