#!/usr/bin/env python3
"""Standalone single-end, both-orientation end-to-end genome mapping."""
import argparse
from genome_common import base_arguments, cli_entry, command, map_run, prepare_run


def build(a,log):
    command(['bowtie2-build','--threads',str(a.threads),a.reference,a.index],log,stdout=log)


def align(a,sam,unmapped,log):
    command(['bowtie2','--very-sensitive','--end-to-end','--seed',str(a.seed),
             '-p',str(a.threads),'-x',a.index,'-U',a.fastq,'--un',str(unmapped),'-S',str(sam)],log)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='operation',required=True)
    base_arguments(sub.add_parser('prepare',help='Build managed index without modifying FASTA'))
    base_arguments(sub.add_parser('map',help='Map trimmed FASTQ using existing index or build from FASTA'),True)
    def run(a):
        if a.operation=='prepare':prepare_run(a,'bowtie2',build)
        else:map_run(a,'bowtie2',align,build)
    cli_entry(p,run)

if __name__=='__main__':main()
