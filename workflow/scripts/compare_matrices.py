#!/usr/bin/env python3
"""Exact keyed comparison of one-based, single-nucleotide integer matrices."""
import argparse
import csv
import json
from pathlib import Path
import re

from bam_count import read_sizes
from core import sha256


def read_matrix(path, sizes, minimum):
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, [])
        if header[:4] != ['chromosome', 'strand', 'start', 'end'] or len(header) < 5:
            raise ValueError('Invalid matrix coordinate header')
        names = header[4:]
        if len(set(names)) != len(names) or any(not name for name in names):
            raise ValueError('Invalid/duplicate matrix sample names')
        sites, totals = {}, [0]*len(names)
        for row in reader:
            if len(row) != len(header) or not all(re.fullmatch(r'\d+', x) for x in [*row[2:4], *row[4:]]):
                raise ValueError('Malformed matrix row or nonnegative-integer violation')
            pos, end = int(row[2]), int(row[3])
            if row[0] not in sizes or row[1] not in ('+', '-') or pos != end or not 1 <= pos <= sizes[row[0]]:
                raise ValueError('Invalid one-based single-nucleotide coordinate')
            key = row[0], row[1], pos
            values = list(map(int, row[4:]))
            if key in sites:
                raise ValueError('Duplicate feature key')
            if sum(values) < minimum:
                raise ValueError('Feature below minimum total-count threshold')
            sites[key] = dict(zip(names, values))
            totals = [a+b for a,b in zip(totals, values)]
    return names, sites, dict(zip(names, totals))


def compare(expected, observed, sizes_path, minimum=0, expected_features=None, expected_sha256=None):
    if type(minimum) is not int or minimum < 0:
        raise ValueError('minimum must be a nonnegative integer')
    if expected_sha256 is not None and sha256(expected) != expected_sha256:
        raise ValueError('Historical matrix SHA-256 mismatch')
    sizes = read_sizes(sizes_path)
    a_names, a, a_totals = read_matrix(expected, sizes, minimum)
    b_names, b, b_totals = read_matrix(observed, sizes, minimum)
    missing, extra = set(a)-set(b), set(b)-set(a)
    mismatches, maximum, examples = 0, 0, []
    by_sample = {name: 0 for name in a_names}
    for key in sorted(set(a)&set(b)):
        for name in a_names:
            if name not in b[key]:
                continue
            delta = b[key][name]-a[key][name]
            if delta:
                mismatches += 1
                by_sample[name] += 1
                maximum = max(maximum, abs(delta))
                if len(examples) < 20:
                    examples.append({'feature': key, 'sample_id': name, 'expected': a[key][name], 'observed': b[key][name]})
    columns_equal = a_names == b_names
    feature_count_ok = expected_features is None or len(a) == len(b) == expected_features
    return {'equal': columns_equal and not missing and not extra and mismatches == 0 and feature_count_ok and a_totals == b_totals,
            'coordinate_definition': 'reference identifier, strand, 1-based position; start=end',
            'sample_order_equal': columns_equal, 'expected_sample_order': a_names, 'observed_sample_order': b_names,
            'expected_features': len(a), 'observed_features': len(b), 'feature_count_requirement_met': feature_count_ok,
            'missing_features': len(missing), 'extra_features': len(extra),
            'missing_examples': sorted(missing)[:20], 'extra_examples': sorted(extra)[:20],
            'unequal_cells': mismatches, 'maximum_absolute_difference': maximum,
            'unequal_cells_by_sample': by_sample, 'cell_examples': examples,
            'expected_sample_totals': a_totals, 'observed_sample_totals': b_totals,
            'sample_totals_equal': a_totals == b_totals, 'integer_tolerance': 0,
            'expected_sha256': sha256(expected), 'observed_sha256': sha256(observed),
            'byte_equal': sha256(expected) == sha256(observed)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected', required=True)
    parser.add_argument('--observed', required=True)
    parser.add_argument('--sizes', required=True, help='Read-only two-column reference sizes or .fai')
    parser.add_argument('--min-total-count', type=int, default=0)
    parser.add_argument('--expected-features', type=int)
    parser.add_argument('--expected-sha256')
    parser.add_argument('--report', required=True, help='New JSON path; never replaces an input or existing result')
    args = parser.parse_args()
    path = Path(args.report).resolve()
    if path in {Path(p).resolve() for p in (args.expected, args.observed, args.sizes)} or path.exists():
        parser.error('Choose a new report path distinct from inputs')
    try:
        result = compare(args.expected, args.observed, args.sizes, args.min_total_count,
                         args.expected_features, args.expected_sha256)
    except (ValueError, OSError) as error:
        result = {'equal': False, 'error': str(error)}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+'\n')
    return 0 if result['equal'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
