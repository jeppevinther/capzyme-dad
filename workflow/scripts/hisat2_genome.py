#!/usr/bin/env python3
"""Standalone single-end splice-aware genome mapping; explicit index preparation."""
import argparse
import json
from pathlib import Path
import subprocess
from genome_common import base_arguments, cli_entry, command, map_run, prepare_run, reference_sizes, index_files


def annotations(path,sizes,kind):
    n=0
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t')
        if kind=='gtf':
            if len(f)!=9:raise ValueError('GTF requires nine columns')
            if f[2]!='exon':continue
            chrom,start,end,strand=f[0],int(f[3])-1,int(f[4])-1,f[6]
        else:
            if len(f)!=4:raise ValueError('Splice/exon TSV requires chrom, start, end, strand')
            chrom,start,end,strand=f[0],int(f[1]),int(f[2]),f[3]
        if chrom not in sizes or strand not in ('+','-') or not 0<=start<=end<sizes[chrom]:
            raise ValueError('Annotation/reference identifiers or coordinates incompatible')
        n+=1
    if not n:raise ValueError('No compatible exon/splice records')


def build(a,log):
    if a.gtf and (a.splice_sites or a.exons):raise ValueError('Use GTF or pre-extracted splice/exon tables')
    out=Path(a.index).parent.parent
    sizes=reference_sizes(a.reference,out/'annotation.sizes.tsv',a.reference_sha256)
    ss,exons=a.splice_sites,a.exons
    if a.gtf:
        annotations(a.gtf,sizes,'gtf')
        ss,exons=str(out/'splice_sites.tsv'),str(out/'exons.tsv')
        for tool,dest in [('hisat2_extract_splice_sites.py',ss),('hisat2_extract_exons.py',exons)]:
            with open(dest,'w') as f:command([tool,a.gtf],log,stdout=f)
        if Path(ss).stat().st_size==0:ss=None  # Valid single-exon GTF has no junctions.
    if exons and not ss:
        raise ValueError('Exon-only/no-junction graph preparation is unsupported; explicitly omit annotations for a linear index or supply compatible splice sites')
    for path,kind in [(ss,'splice'),(exons,'exon')]:
        if path:annotations(path,sizes,kind)
    args=['hisat2-build','-p',str(a.threads)]
    if ss:args+=['--ss',ss]
    if exons:args+=['--exon',exons]
    command(args+[a.reference,a.index],log,stdout=log)


def align(a,sam,unmapped,log):
    # HISAT2's Perl wrapper shell-quotes some paths incompletely. Give it only
    # simple relative aliases in an output-local working directory, preserving
    # the wrapper's exact --un semantics without editing vendor code.
    work=sam.parent
    for source in index_files(a.index,'hisat2'):
        suffix=str(source)[len(str(Path(a.index).resolve())):]
        (work/('genome'+suffix)).symlink_to(source)
    reads='reads.fq.gz' if a.fastq.endswith('.gz') else 'reads.fq'
    (work/reads).symlink_to(Path(a.fastq).resolve())
    log.write(json.dumps({'working_directory':str(work),'readonly_index_alias':str(Path(a.index).resolve()),
                          'readonly_fastq_alias':str(Path(a.fastq).resolve())})+'\n');log.flush()
    command(['hisat2','--very-sensitive','--no-softclip','--seed',str(a.seed),
             '-p',str(a.threads),'-x','genome','-U',reads,'--un','unmapped.fastq','-S',sam.name],log,cwd=work)
    (work/'unmapped.fastq').replace(unmapped)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='operation',required=True)
    prep=sub.add_parser('prepare',help='Explicit small/large genome index preparation; annotations optional')
    base_arguments(prep);prep.add_argument('--gtf');prep.add_argument('--splice-sites');prep.add_argument('--exons')
    base_arguments(sub.add_parser('map',help='Map trimmed FASTQ using a pre-existing index'),True)
    def run(a):
        if a.operation=='prepare':prepare_run(a,'hisat2',build)
        elif not a.index:raise ValueError('HISAT2 map requires --index; use prepare separately')
        else:map_run(a,'hisat2',align)
    cli_entry(p,run)

if __name__=='__main__':main()
