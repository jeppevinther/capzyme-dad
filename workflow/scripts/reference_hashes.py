"""Observed SHA-256 provenance and optional, binding expected identities."""
import json
from pathlib import Path
import re
from core import sha256


def expected_digest(value, label):
    if value is None: return None
    if not isinstance(value,str) or not re.fullmatch('[0-9a-fA-F]{64}',value):
        raise ValueError(label+' must be null or a 64-character hexadecimal SHA-256')
    return value.lower()


def observed_file(path, expected=None, label='Reference SHA-256'):
    pin=expected_digest(expected,label)
    if not Path(path).is_file(): raise ValueError('Missing reference/index input: '+str(path))
    observed=sha256(path)
    if pin is not None and pin!=observed: raise ValueError(label+' mismatch: '+str(path))
    return observed


def check_manifest(path, reference, prefix, kind):
    """Trusted manifests are complete earlier records, never partial expectations."""
    from genome_common import index_files
    try: manifest=json.loads(Path(path).read_text())
    except (json.JSONDecodeError,UnicodeError) as exc: raise ValueError('Malformed trusted manifest') from exc
    if not isinstance(manifest,dict): raise ValueError('Trusted manifest must be an object')
    files=index_files(prefix,kind)
    observed={str(p)[len(str(Path(prefix).resolve())):]:sha256(p) for p in files}
    parts=manifest.get('parts')
    if not isinstance(parts,dict) or set(parts)!=set(observed): raise ValueError('Trusted manifest must contain the complete index family')
    normalized={k:expected_digest(v,'Trusted index hash') for k,v in parts.items()}
    fasta=expected_digest(manifest.get('reference_sha256'),'Trusted FASTA hash')
    if fasta is None or any(v is None for v in normalized.values()): raise ValueError('Trusted manifest hashes cannot be null or omitted')
    if type(manifest.get('schema_version')) is not int or manifest.get('schema_version')!=1 or manifest.get('aligner')!=kind or fasta!=sha256(reference) or normalized!=observed:
        raise ValueError('Reference/index manifest checksum or aligner mismatch')
    return manifest
