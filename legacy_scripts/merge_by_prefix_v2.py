#!/usr/bin/env python3
"""
Merge Ensembl ncRNA FASTA records by identity of the first K bases.

Default behavior:
- Groups sequences with identical first K nucleotides, case-insensitive.
  If a sequence is shorter than K, its full length is used as the key.
- Drops records whose biotype is listed in --drop-biotypes.
- Drops records whose description field mentions the word "pseudogene".
- Representative sequence for each group = longest sequence in the group.
  Ties are resolved by taking the first sequence encountered in the input FASTA.
- Merged FASTA IDs are short family-style names:
      <family_or_symbol>__n<collapsed_member_count>
  Example:
      RNU2__n37
      SNORD14C__n3
      ENST00000612345.1__n1
- If two groups would receive the same short ID, the later groups get a short
  numeric suffix:
      RNU2__n37
      RNU2__n37__g2
      RNU2__n37__g3
  This preserves unique FASTA record names without using long member lists.
- Full collapse information is written to the membership TSV.

Example:
  ./merge_by_prefix_v2.py Homo_sapiens.GRCh38.ncrna.fa \
      -o merged_ncRNA_70bp \
      -k 70 \
      --drop-biotypes retained_intron,pseudogene,snRNA_pseudogene,snoRNA_pseudogene,miRNA_pseudogene
"""

from __future__ import annotations
import argparse
import gzip
import io
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

FASTA_WRAP = 60

# -------- I/O helpers --------

def open_maybe_gzip(path: str | Path) -> io.TextIOBase:
    """Open plain-text or gzip-compressed input in text mode."""
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    return open(p, "r", encoding="utf-8")


def _cleanup_seq(seq: str) -> str:
    """Remove whitespace and uppercase the sequence for case-insensitive matching."""
    return re.sub(r"\s+", "", seq).upper()


def parse_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (header_without_>, cleaned_sequence) from a FASTA file."""
    with open_maybe_gzip(path) as fh:
        header = None
        chunks: List[str] = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, _cleanup_seq("".join(chunks))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, _cleanup_seq("".join(chunks))


# -------- header parsing --------

GENE_SYMBOL_PATTERNS = [
    r'(?:^|\s)gene_symbol[:=]\s*"?([A-Za-z0-9._-]+)"?',
    r'(?:^|\s)gene_symbol\s+"([^"]+)"',
]

# Anchored to field boundaries so "biotype" does not also match inside
# "gene_biotype" and "transcript_biotype".
BIOTYPE_PATTERNS = [
    r'(?:^|\s)gene_biotype[:=]\s*"?([A-Za-z0-9_+-]+)"?',
    r'(?:^|\s)transcript_biotype[:=]\s*"?([A-Za-z0-9_+-]+)"?',
    r'(?:^|\s)biotype[:=]\s*"?([A-Za-z0-9_+-]+)"?',
    r'(?:^|\s)gene_biotype\s+"([^"]+)"',
    r'(?:^|\s)transcript_biotype\s+"([^"]+)"',
    r'(?:^|\s)biotype\s+"([^"]+)"',
]

# Look specifically in the description field and match whole-word pseudogene.
PSEUDO_DESC_PAT = re.compile(r'\bdescription\b[^>]*?\bpseudogene\b', re.IGNORECASE)


def extract_ids(header: str) -> Tuple[str, str | None]:
    """Return primary FASTA ID and optional gene_symbol from an Ensembl-style header."""
    primary_id = header.split()[0]
    gene_symbol = None
    for pat in GENE_SYMBOL_PATTERNS:
        m = re.search(pat, header)
        if m:
            gene_symbol = m.group(1)
            break
    return primary_id, gene_symbol


def extract_biotypes(header: str) -> List[str]:
    """Return unique biotypes from an Ensembl-style header, preserving discovery order."""
    bios: List[str] = []
    seen = set()
    for pat in BIOTYPE_PATTERNS:
        for m in re.finditer(pat, header):
            b = m.group(1)
            if b not in seen:
                bios.append(b)
                seen.add(b)
    return bios or ["(none)"]


def is_pseudogene_description(header: str) -> bool:
    """Return True if the header's description field mentions 'pseudogene'."""
    return bool(PSEUDO_DESC_PAT.search(header))


# -------- merged ID helpers --------

def simplify_symbol(symbol: str) -> str:
    """
    Convert common small-RNA copy symbols to family-style labels.

    Examples:
      RNU2-1   -> RNU2
      RNU6-74  -> RNU6
      RNU4-2   -> RNU4

    Other symbols are returned unchanged.
    """
    # U-RNA/snRNA copy naming commonly uses RNU<number>-<copy>.
    m = re.match(r'^(RNU\d+[A-Za-z]?)-\d+[A-Za-z]*$', symbol)
    if m:
        return m.group(1)
    return symbol


def clean_id_component(text: str, max_len: int = 60) -> str:
    """Make a FASTA-ID-safe, compact component."""
    cleaned = re.sub(r'[^A-Za-z0-9_.-]+', '_', text).strip('_')
    return (cleaned[:max_len].strip('_') or "ncRNA")


def make_base_merged_id(prefix: str, member_count: int) -> str:
    """Create the requested short merged ID before collision handling."""
    return f"{clean_id_component(prefix)}__n{member_count}"


def uniquify_id(base_id: str, used: Dict[str, int]) -> str:
    """
    Ensure FASTA IDs are unique.

    The first occurrence keeps base_id. Later collisions get __g2, __g3, ...
    This avoids long hashes/member lists while preserving valid unique reference names.
    """
    used[base_id] += 1
    if used[base_id] == 1:
        return base_id
    return f"{base_id}__g{used[base_id]}"


# -------- merging --------

def wrap_fasta(seq: str, width: int = FASTA_WRAP) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def merge_by_prefix(
    fasta_path: str | Path,
    k: int = 70,
    drop_biotypes: set[str] | None = None,
) -> Tuple[List[Tuple[str, str]], List[Dict[str, str]], Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Return merged FASTA records, membership rows, pre-filter biotype counts,
    dropped-reason counts, technical skip counts, and record-level summary stats.
    """
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    biotype_counts_pre: Counter = Counter()
    biotype_counts_dropped: Counter = Counter()
    skip_counts: Counter = Counter()
    stats: Counter = Counter()
    drop_biotypes = {b.lower() for b in (drop_biotypes or set())}

    for header, seq in parse_fasta(fasta_path):
        stats['records_seen'] += 1
        if not seq:
            skip_counts["empty_sequence"] += 1
            stats['records_skipped'] += 1
            continue

        bios = extract_biotypes(header)
        tid, symbol = extract_ids(header)

        for b in bios:
            biotype_counts_pre[b] += 1

        drop_due_to_named = any(b.lower() in drop_biotypes for b in bios)
        drop_due_to_pseudogene_desc = is_pseudogene_description(header)

        if drop_due_to_named or drop_due_to_pseudogene_desc:
            if drop_due_to_named:
                for b in bios:
                    if b.lower() in drop_biotypes:
                        biotype_counts_dropped[b] += 1
            if drop_due_to_pseudogene_desc:
                biotype_counts_dropped["__pseudogene_description__"] += 1
            stats['records_dropped'] += 1
            continue

        key = seq[:k]
        stats['records_kept'] += 1
        groups[key].append({
            "header": header,
            "seq": seq,
            "len": len(seq),
            "id": tid,
            "symbol": symbol or "",
            "family_symbol": simplify_symbol(symbol) if symbol else "",
            "biotypes": ",".join(bios),
        })

    merged_records: List[Tuple[str, str]] = []
    membership_rows: List[Dict[str, str]] = []
    used_merged_ids: Dict[str, int] = defaultdict(int)

    for key, members in groups.items():
        ids = [str(m["id"]) for m in members]
        symbols = [str(m["symbol"]) for m in members if m["symbol"]]
        family_symbols = [str(m["family_symbol"]) for m in members if m["family_symbol"]]

        # Choose family-style merged ID prefix.
        if family_symbols:
            prefix = Counter(family_symbols).most_common(1)[0][0]
        else:
            prefix = ids[0]

        base_id = make_base_merged_id(prefix, len(members))
        merged_id = uniquify_id(base_id, used_merged_ids)

        # Representative = longest; tie -> first encountered.
        _, rep = max(enumerate(members), key=lambda t: (int(t[1]["len"]), -t[0]))
        rep_seq = str(rep["seq"])

        merged_records.append((merged_id, rep_seq))
        membership_rows.append({
            "merged_id": merged_id,
            "base_merged_id": base_id,
            "key_firstK": key,
            "member_count": str(len(members)),
            "representative_id": str(rep["id"]),
            "representative_len": str(rep["len"]),
            "representative_symbol": str(rep["symbol"]),
            "family_symbol": prefix,
            "member_ids": ",".join(ids),
            "member_symbols": ",".join(sorted(set(symbols))) if symbols else "",
            "member_family_symbols": ",".join(sorted(set(family_symbols))) if family_symbols else "",
        })

    return (
        merged_records,
        membership_rows,
        dict(biotype_counts_pre),
        dict(biotype_counts_dropped),
        dict(skip_counts),
        dict(stats),
    )


# -------- writers & CLI --------

def write_fasta(records: Iterable[Tuple[str, str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as fo:
        for rid, seq in records:
            fo.write(f">{rid}\n{wrap_fasta(seq)}\n")


def write_tsv(rows: Iterable[Dict[str, str]], path: str | Path) -> None:
    cols = [
        "merged_id",
        "base_merged_id",
        "key_firstK",
        "member_count",
        "representative_id",
        "representative_len",
        "representative_symbol",
        "family_symbol",
        "member_ids",
        "member_symbols",
        "member_family_symbols",
    ]
    with open(path, "w", encoding="utf-8") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(r.get(c, "") for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Merge ncRNA FASTA by first-K bases; use short family-style IDs; drop pseudogenes by description and optional biotypes."
    )
    ap.add_argument("fasta", help="Input Ensembl ncRNA FASTA (.fa/.fa.gz).")
    ap.add_argument("-o", "--out-prefix", default="merged_ncRNA_70bp", help="Output prefix.")
    ap.add_argument("-k", "--key-length", type=int, default=70, help="Prefix length for grouping (default 70).")
    ap.add_argument(
        "--drop-biotypes",
        type=str,
        default="",
        help="Comma-separated list of biotypes to drop (case-insensitive), e.g. 'retained_intron,pseudogene'.",
    )
    # Kept for backwards compatibility with older commands. It no longer affects IDs,
    # because member lists are now stored only in the TSV.
    ap.add_argument(
        "--id-suffix-count",
        type=int,
        default=0,
        help=argparse.SUPPRESS,
    )
    args = ap.parse_args()

    drop_set = {b.strip() for b in args.drop_biotypes.split(",") if b.strip()}

    merged, rows, pre, dropped, skipped, stats = merge_by_prefix(
        args.fasta,
        k=args.key_length,
        drop_biotypes=drop_set,
    )

    out_fa = f"{args.out_prefix}.fa"
    out_tsv = f"{args.out_prefix}.tsv"
    write_fasta(merged, out_fa)
    write_tsv(rows, out_tsv)

    total_records_seen = stats.get('records_seen', 0)
    total_dropped_records = stats.get('records_dropped', 0)
    total_dropped_reasons = sum(dropped.values())

    print("\n[Biotype counts in headers] (top 20)")
    for bt, n in sorted(pre.items(), key=lambda x: -x[1])[:20]:
        print(f"  {bt:20s}  {n:,}")

    if total_dropped_records:
        print("\n[Dropped records]")
        for bt, n in sorted(dropped.items(), key=lambda x: -x[1]):
            label = "pseudogene(desc)" if bt == "__pseudogene_description__" else bt
            print(f"  {label:20s}  {n:,}")
        print(f"  Total dropped records: {total_dropped_records:,} / {total_records_seen:,}")
        if total_dropped_reasons != total_dropped_records:
            print(f"  Note: dropped-reason count is {total_dropped_reasons:,} because records can match >1 reason.")
    else:
        print("\n[Dropped records] None matched criteria.")

    if skipped:
        print("\n[Skipped records (technical)]")
        for why, n in sorted(skipped.items(), key=lambda x: -x[1]):
            print(f"  {why:20s}  {n:,}")

    print("\n[Record summary]")
    print(f"  Records seen:    {stats.get('records_seen', 0):,}")
    print(f"  Records kept:    {stats.get('records_kept', 0):,}")
    print(f"  Records dropped: {stats.get('records_dropped', 0):,}")

    print("\n[Merge]")
    print(f"  Groups (unique by first {args.key_length}): {len(rows):,}")
    singles = sum(1 for r in rows if r["member_count"] == "1")
    print(f"  Groups with single member: {singles:,}")

    collisions = sum(1 for r in rows if r["merged_id"] != r["base_merged_id"])
    if collisions:
        print(f"  Short-ID collisions resolved with __gN suffixes: {collisions:,}")

    print(f"\n[Outputs]\n  FASTA: {out_fa}\n  TSV:   {out_tsv}")


if __name__ == "__main__":
    main()
