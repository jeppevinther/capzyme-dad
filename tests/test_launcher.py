"""Exercise the beginner CLI in a fresh external working directory."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT/'capzyme.py'
TOOLS = all(shutil.which(t) for t in ('snakemake','cutadapt','bowtie2','bowtie2-build',
                                     'bowtie2-inspect','samtools','bedtools'))


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='capzyme beginner ')
        self.addCleanup(self.temp.cleanup)
        self.p = Path(self.temp.name).resolve()

    def call(self, *args, ok=True, env=None):
        r = subprocess.run([sys.executable,str(CLI),*map(str,args)],cwd=self.p,
                           text=True,capture_output=True,env=env)
        if ok:self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        else:self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
        return r

    def test_example_validation_paths_and_nonoverwriting(self):
        demo = self.p/'demo'
        self.call('example','--directory',demo)
        cfg = demo/'config.yaml'
        result = self.call('check','--config',cfg)
        self.assertIn('B, A, empty',result.stdout)
        self.assertIn(str(demo),result.stdout)
        self.assertFalse((demo/'run').exists())
        before = cfg.read_bytes()
        self.call('example','--directory',demo,ok=False)
        self.assertEqual(cfg.read_bytes(),before)
        # Relative paths bind to the launcher's external workdir, not the caller.
        data=json.loads(cfg.read_text());data['samples']='samples.tsv'
        data['reference']['ncrna_fasta']='ncrna.fa'
        data['reference']['viral_fastas']=['viral.fa'];data['output_dir']='run'
        cfg.write_text(json.dumps(data))
        self.call('check','--config',cfg)
        self.assertIn('does not exist',self.call('check','--config',cfg,'--workdir',self.p,ok=False).stderr)
        self.assertIn('outside',self.call('check','--config',cfg,'--workdir',ROOT,ok=False).stderr)
        self.call('example','--directory',ROOT/'unwanted-tutorial',ok=False)
        self.assertFalse((ROOT/'unwanted-tutorial').exists())

    def test_invalid_yaml_route_dependency_and_resource_messages(self):
        self.call('example','--directory',self.p/'demo')
        cfg=self.p/'demo/config.yaml'
        self.call('check','--config',cfg,'--route','sequential',ok=False)
        self.call('preview','--config',cfg,'--cores','0',ok=False)
        cfg.write_text('schema_version: [unclosed')
        self.assertIn('Invalid YAML',self.call('check','--config',cfg,ok=False).stderr)
        import os
        empty_env=dict(os.environ,PATH=str(self.p/'missing-tools'))
        result=self.call('doctor',env=empty_env,ok=False)
        self.assertIn('MISSING snakemake',result.stdout)
        self.assertIn('MISSING bowtie2',result.stdout)

    @unittest.skipUnless(TOOLS,'Full RNA workflow tools required')
    def test_documented_preview_run_exact_results_and_noop_restart(self):
        self.call('doctor')
        demo=self.p/'example';self.call('example','--directory',demo)
        cfg=demo/'config.yaml'
        self.call('check','--config',cfg)
        self.call('preview','--config',cfg,'--cores','4','--memory-mb','16384')
        self.assertFalse((demo/'run/counts/5ends.raw.tsv').exists())
        self.call('run','--config',cfg,'--cores','4','--memory-mb','16384')
        matrix=demo/'run/counts/5ends.raw.tsv'
        self.assertEqual(matrix.read_bytes(),(demo/'expected.tsv').read_bytes())
        self.assertEqual(len(matrix.read_text().splitlines()),5)
        filtered=demo/'run/counts/5ends.min_total_10.tsv'
        self.assertEqual(filtered.read_text().splitlines()[1],'RNU2__n2\t+\t1\t1\t12\t0\t0')
        import hashlib
        outputs={p:(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns)
                 for p in (demo/'run').rglob('*') if p.is_file()}
        again=self.call('run','--config',cfg,'--cores','4','--memory-mb','16384')
        self.assertIn('Nothing to be done',again.stdout+again.stderr)
        for path,expected in outputs.items():
            self.assertEqual((hashlib.sha256(path.read_bytes()).hexdigest(),path.stat().st_mtime_ns),expected)
        # Adding tracks to compatible completed RNA counts must not rerun upstream.
        data=json.loads(cfg.read_text())
        data['bigwigs']={'rna':{'methods':{'raw':{},'rpm':{}},'cohort':['B','A'],'strands':'plus'}}
        cfg.write_text(json.dumps(data))
        preview=self.call('preview','--config',cfg)
        for name in ['trim_fastq','map_filter_sort','count_5ends']:
            self.assertNotIn('rule '+name+':',preview.stdout+preview.stderr)
        self.call('run','--config',cfg)
        for path,expected in outputs.items():
            if 'report' in path.relative_to(demo/'run').parts or path.name=='report.log':
                continue  # Final report changes with track requests; scientific outputs stay fixed.
            self.assertEqual((hashlib.sha256(path.read_bytes()).hexdigest(),path.stat().st_mtime_ns),expected)
        report=json.loads((demo/'run/report/report_status.json').read_text())
        self.assertEqual(report['status'],'complete_from_records')
        # Deleting the report schedules reporting only, with no dependency cycle.
        (demo/'run/report/report.html').unlink()
        dry=self.call('preview','--config',cfg)
        self.assertIn('rule run_report:',dry.stdout+dry.stderr)
        for rule in ['trim_fastq','map_filter_sort','count_5ends','bigwig_rna_raw','bigwig_rna_rpm']:
            self.assertNotIn('rule '+rule+':',dry.stdout+dry.stderr)
        self.call('run','--config',cfg)
        self.assertEqual(json.loads((demo/'run/report/analysis_context.json').read_text())['route'],'rna')
        self.assertTrue((demo/'run/bigwigs/rna/rpm/tracks/A.plus.bw').is_file())
        # The adapter's other routes must reach the existing validated entry points.
        import hashlib
        reference=demo/'run/reference/combined.fa'
        digest=hashlib.sha256(reference.read_bytes()).hexdigest()
        manifest=demo/'bams.tsv'
        manifest.write_text('sample_id\tbam\n'+''.join(f'{n}\t{demo/"run/bam"/(n+".accepted.bam")}\n' for n in ['B','A','empty']))
        bam_config={'schema_version':1,'input_mode':'existing_bam','samples':str(manifest),
                    'output_dir':str(demo/'bam-counts'),
                    'reference':{'fasta':str(reference),'fasta_sha256':digest},'counting':{'strand':'both'},'bigwigs':{'rna':{'methods':{'raw':{},'rpm':{}},'cohort':['B','A']}}}
        bam_path=demo/'bam.json';bam_path.write_text(json.dumps(bam_config))
        for operation in ['check','preview','run']:
            self.call(operation,'--route','bam','--config',bam_path)
        self.assertEqual((demo/'bam-counts/counts/5ends.raw.tsv').read_bytes(),matrix.read_bytes())
        reads=demo/'genome-samples.tsv'
        reads.write_text(f'sample_id\tfastq\nB\t{demo/"run/trimmed/B.fastq.gz"}\n')
        genome={'aligner':'bowtie2','reference':str(reference),'reference_sha256':digest,
                'index_prefix':str(demo/'run/reference/bowtie2/index'),'no_decoys':True,'threads':4,'mem_mb':4000}
        genome_path=demo/'genome.json'
        genome_path.write_text(json.dumps({'schema_version':1,'samples':str(reads),
                                          'output_dir':str(demo/'genome-only'),'genome':genome,'bigwigs':{'genome':{'methods':{'raw':{},'rpm':{}}}}}))
        for operation in ['check','preview','run']:
            self.call(operation,'--route','genome','--config',genome_path)
        self.assertEqual((demo/'genome-only/counts/genome.5ends.raw.tsv').read_text(),
                         'chromosome\tstrand\tstart\tend\tB\nRNU2__n2\t+\t1\t1\t12\n')
        sequential=json.loads(cfg.read_text());sequential['genome']=genome
        sequential['bigwigs']['genome']={'methods':{'raw':{}},'strands':'both'}
        sequential['output_dir']=str(demo/'sequential')
        sequential_path=demo/'sequential.json';sequential_path.write_text(json.dumps(sequential))
        for operation in ['check','preview','run']:
            self.call(operation,'--route','sequential','--config',sequential_path)
        report=json.loads((demo/'sequential/qc/B.sequential.json').read_text())
        self.assertTrue((demo/'sequential/bigwigs/genome/raw/tracks/A.minus.bw').is_file())
        self.assertTrue(report['handoff_sequence_quality_verified'])
        self.assertEqual(report['accepted_read_id_overlap'],0)


if __name__ == '__main__':unittest.main()
