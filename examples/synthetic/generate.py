"""Generate a deterministic artificial RNA-only tutorial; no research data."""
import json
from pathlib import Path
import random

ADAPTER = 'AGATCGGAAGAGCACACGTCT'


def sequence(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def generate(directory):
    directory = Path(directory).resolve()
    repo = Path(__file__).resolve().parents[2]
    if directory == repo or repo in directory.parents:
        raise ValueError('Synthetic example must be outside the source checkout')
    if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
        raise ValueError('Synthetic example directory must be new or empty')
    directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(713)
    rna, virus = sequence(rng, 160), sequence(rng, 140)
    minus = virus.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    (directory/'ncrna.fa').write_text(f'>short gene_symbol:RNU2-1\n{rna[:100]}\n>long gene_symbol:RNU2-2\n{rna}\n')
    (directory/'viral.fa').write_text(f'>virus_plus\n{virus}\n>virus_minus\n{minus}\n')
    reads = {'B': [rna[:45]]*12,
             'A': [rna[10:55]]*3+[virus[9:54], minus[20:65]], 'empty': ['']}
    for sample, sequences in reads.items():
        with (directory/f'{sample}.fastq').open('w') as fh:
            for i, seq in enumerate(sequences):
                full = seq+ADAPTER
                fh.write(f'@{sample}_{i}\n{full}\n+\n'+'I'*len(full)+'\n')
    (directory/'samples.tsv').write_text('sample_id\tfastq\n'+''.join(f'{name}\t{directory/(name+".fastq")}\n' for name in reads))
    # JSON is a YAML subset. All data/output paths are absolute and external.
    config = {'schema_version': 1, 'samples': str(directory/'samples.tsv'),
              'output_dir': str(directory/'run'),
              'reference': {'ncrna_fasta': str(directory/'ncrna.fa'),
                            'viral_fastas': [str(directory/'viral.fa')], 'prefix_length': 30, 'drop_biotypes': []},
              'preprocessing': {'adapter': ADAPTER, 'quality': 20, 'minimum_length': 15, 'nextseq_trim': 20},
              'mapping': {'norc': True, 'preset': 'very-sensitive', 'seed': 0, 'threads': 4},
              'counting': {'mapq': 20, 'strand': 'both', 'exclude_duplicates': False, 'min_total_count': 10}}
    path = directory/'config.yaml'
    path.write_text(json.dumps(config, indent=2)+'\n')
    (directory/'expected.tsv').write_text(
        'chromosome\tstrand\tstart\tend\tB\tA\tempty\n'
        'RNU2__n2\t+\t1\t1\t12\t0\t0\n'
        'RNU2__n2\t+\t11\t11\t0\t3\t0\n'
        'virus_minus\t+\t21\t21\t0\t1\t0\n'
        'virus_plus\t+\t10\t10\t0\t1\t0\n')
    return path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    print(generate(parser.parse_args().directory))
