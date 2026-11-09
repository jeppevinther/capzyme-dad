"""Read-only validation of existing FASTA and Bowtie2 indexes."""
import json
from pathlib import Path
import re
import subprocess
import tempfile
from core import records, sha256, validate_fasta
from reference_hashes import expected_digest

PARTS = ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2')


def fixed_inputs(reference):
    fmt = reference.get('index_format', 'bt2')
    if fmt not in ('bt2', 'bt2l'):
        raise ValueError('reference.index_format must be bt2 or bt2l')
    prefix = reference['index_prefix']
    suffixes = [part + '.' + fmt for part in PARTS]
    return suffixes, [prefix + suffix for suffix in suffixes]


def check_fixed_hashes(reference):
    suffixes, files = fixed_inputs(reference)
    pins = reference.get('index_sha256')
    if pins is None: pins = {}
    if not isinstance(pins,dict) or set(pins)-set(suffixes):
        raise ValueError('reference.index_sha256 must map known index suffixes to expected hashes')
    # Validate the entire supplied specification before reading any artifact.
    normalized={key:expected_digest(value,'Expected index SHA-256') for key,value in pins.items()}
    if any(value is None for value in normalized.values()):
        raise ValueError('Individual expected index hashes cannot be null; omit unpinned suffixes')
    fasta_pin=expected_digest(reference.get('fasta_sha256'),'Expected FASTA SHA-256')
    expected = [(reference['fasta'], fasta_pin)]
    expected += [(path, normalized.get(suffix)) for suffix, path in zip(suffixes, files)]
    other = 'bt2l' if reference.get('index_format', 'bt2') == 'bt2' else 'bt2'
    if any(Path(reference['index_prefix'] + part + '.' + other).exists() for part in PARTS):
        raise ValueError('Ambiguous index prefix: both bt2 and bt2l families are present')
    observed = {}
    for path, pin in expected:
        if not Path(path).is_file():
            raise ValueError('Missing fixed-reference input: ' + str(path))
        digest = sha256(path)
        if pin is not None and digest != pin:
            raise ValueError('SHA-256 mismatch for fixed-reference input: ' + str(path))
        observed[str(path)] = digest
    return observed


def verify_fixed(reference, sizes, report):
    before = check_fixed_hashes(reference)
    validate_fasta(reference['fasta'])
    expected = list(records(reference['fasta']))
    # Inspect into an output-local temporary file, never into the input directory.
    with tempfile.NamedTemporaryFile(dir=Path(report).parent, suffix='.inspect.tmp') as tmp:
        proc = subprocess.run(['bowtie2-inspect', reference['index_prefix']],
                              stdout=tmp, stderr=subprocess.PIPE, text=True)
        if proc.returncode:
            raise ValueError('bowtie2-inspect failed: ' + proc.stderr)
        tmp.flush()
        actual = list(records(tmp.name))
    # Bowtie2 stores ambiguous IUPAC bases as N in reconstructed sequences.
    masked = [(name, re.sub('[^ACGT]', 'N', sequence)) for name, sequence in expected]
    if actual != masked:
        raise ValueError('Existing Bowtie2 index does not match the fixed FASTA identifiers/order/sequences')
    after = check_fixed_hashes(reference)
    if before != after:
        raise ValueError('Fixed-reference inputs changed during verification')
    with open(sizes, 'w') as fh:
        for name, sequence in expected:
            fh.write(f'{name}\t{len(sequence)}\n')
    version = subprocess.check_output(['bowtie2-inspect', '--version'], text=True).splitlines()[0]
    result = {'mode': 'fixed', 'inputs_sha256': after, 'sequence_count': len(expected),
              'total_bases': sum(len(seq) for _, seq in expected),
              'index_format': reference.get('index_format', 'bt2'),
              'index_matches_fasta': True,
              'expected_hashes': {'fasta_sha256': reference.get('fasta_sha256'), 'index_sha256': reference.get('index_sha256')},
              'hash_trust': 'observed hashes; supplied expectations enforced where present; unpinned observations are not historical identity proof', 'bowtie2_inspect_version': version,
              'ambiguity_comparison': 'non-ACGT reference symbols masked to N for index comparison only'}
    Path(report).write_text(json.dumps(result, indent=2) + '\n')
    return result
