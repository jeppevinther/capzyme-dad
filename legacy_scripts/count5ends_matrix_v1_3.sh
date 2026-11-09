#!/usr/bin/env bash
set -euo pipefail
# 2026-06-04 (v1.3) count5ends_matrix
# Jeppe Vinther (+ChatGPT)
# =================================================================
# Purpose
# -------
# Build a DESeq2-ready matrix of 5' end counts (single-end libraries).
# One row per genomic position (per strand), with a count column per sample.
#
# For CapZyme-seq analysis, the informative position is the
# 5′ base of each aligned read. We count those per base, strand-specific.
#
# Output format (TSV)
# -------------------
# chromosome  strand  start  end  Sample1  Sample2  ...
# • Coordinates are 1-based (featureCounts-like). We set end == start.
# • “strand” is '+' or '-'.
#
# Filtering policy
# ----------------
# • Strict primary mapped alignments only:
#     samtools -F 2308
#   which excludes:
#     0x4   unmapped
#     0x100 secondary
#     0x800 supplementary
#
# • MAPQ filter:
#     default: --mapq 20
#   Change with --mapq N, or use --mapq 0 to disable MAPQ filtering.
#
# • Duplicates (0x400) are kept by default.
#   If you have UMIs, deduplicate first with UMI-tools, then run this script.
#
# Temporary files
# ---------------
# • By default, temporary files are created in a per-run directory under
#   the current working directory:
#
#       ./5end_tmp.XXXXXXXX/
#
# • Use --tmp-dir DIR to place that temporary directory somewhere else,
#   e.g. a project scratch folder.
# • The script exports TMPDIR and passes samtools sort -T so that common
#   temporary files also stay inside this directory rather than /tmp.
# • Temporary files are deleted on normal exit, Ctrl-C, and TERM.
# • Set KEEP_TMP=1 to keep the temp directory for debugging.
#
# Dependencies
# ------------
# samtools, bedtools (in $PATH)
#
# Usage
# -----
#   count5ends_matrix.sh -g reference.fa \
#     [--min-count N] [--mapq N] [--strand plus|minus|both] [--threads N] \
#     [--tmp-dir DIR] [-o out.tsv] sample1.bam [sample2.bam ...]
#
# Examples
# --------
#   # Default: MAPQ >= 20, strict primary only
#   count5ends_matrix.sh -g reference.fa -o counts.tsv sample1.bam sample2.bam
#
#   # Disable MAPQ filtering
#   count5ends_matrix.sh -g reference.fa --mapq 0 -o counts.tsv sample1.bam sample2.bam
#
#   # More stringent MAPQ
#   count5ends_matrix.sh -g reference.fa --mapq 30 -o counts.tsv sample1.bam sample2.bam
#
# Notes
# -----
# • Use --strand plus when reads were mapped to plus-only RNA references
#   (e.g. transcriptome-only mapping where minus strand has no biological signal).
# • The script merges samples using bedtools unionbedg; if only one sample has
#   signal on a strand, it bypasses unionbedg and still emits a full matrix
#   (filling 0 in other sample columns).
# =================================================================

# ---------- Defaults / CLI options ----------
min_count=0          # row filter: keep only rows with total >= this across samples
mapq=20              # alignment filter: keep alignments with MAPQ >= this
output=""            # output path; if empty → stdout
genome=""            # reference FASTA (to read chrom sizes via .fai)
strand_choice="both" # 'plus', 'minus', or 'both'
threads=1            # samtools sort/index/view threads
tmp_base="."         # parent directory for per-run temporary directory

# Simple option parser
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genome)   genome="$2"; shift 2 ;;
    --min-count)   min_count="$2"; shift 2 ;;
    --mapq)        mapq="$2"; shift 2 ;;
    --strand)      strand_choice="$2"; shift 2 ;;
    --threads)     threads="$2"; shift 2 ;;
    --tmp-dir)     tmp_base="$2"; shift 2 ;;
    -o|--output)   output="$2"; shift 2 ;;
    -h|--help)     sed -n '2,260p' "$0"; exit 0 ;;
    --)            shift; break ;;
    -*)            printf 'Error: unknown option %s\n' "$1" >&2; exit 1 ;;
    *)             break ;;
  esac
done

# Remaining args are BAM files
if [[ -z "${genome:-}" || $# -lt 1 ]]; then
  printf 'Error: provide -g <reference.fa> and at least one BAM\n' >&2
  exit 1
fi
bam_files=("$@")

# Validate option values
case "$strand_choice" in plus|minus|both) ;; *)
  printf 'Error: --strand must be one of: plus | minus | both\n' >&2
  exit 1
  ;;
esac

[[ "$threads" =~ ^[0-9]+$ ]] && [[ "$threads" -ge 1 ]] || {
  printf 'Error: --threads must be >= 1\n' >&2
  exit 1
}

[[ "$min_count" =~ ^[0-9]+$ ]] || {
  printf 'Error: --min-count must be a non-negative integer\n' >&2
  exit 1
}

[[ "$mapq" =~ ^[0-9]+$ ]] || {
  printf 'Error: --mapq must be a non-negative integer\n' >&2
  exit 1
}

# ---------- Tool availability ----------
command -v samtools >/dev/null || {
  printf 'Error: samtools not found\n' >&2
  exit 127
}

command -v bedtools >/dev/null || {
  printf 'Error: bedtools not found\n' >&2
  exit 127
}

# ---------- FASTA index → chromosome sizes ----------
# We read the .fai FASTA index to get contig lengths used later by bedtools unionbedg.
[[ -f "$genome" ]] || {
  printf 'Error: FASTA not found: %s\n' "$genome" >&2
  exit 1
}

[[ -f "${genome}.fai" ]] || samtools faidx "$genome" 1>&2

# ---------- Temporary directory ----------
# Create one per-run temp directory under --tmp-dir, default current directory.
# This avoids filling the server-wide /tmp and makes cleanup simple.
mkdir -p "$tmp_base" || {
  printf 'Error: could not create/access --tmp-dir: %s\n' "$tmp_base" >&2
  exit 1
}

tmpdir=$(mktemp -d "${tmp_base%/}/5end_tmp.XXXXXXXX") || {
  printf 'Error: failed to create temporary directory under %s\n' "$tmp_base" >&2
  exit 1
}

# Many tools respect TMPDIR; samtools sort is also given -T explicitly below.
export TMPDIR="$tmpdir"

cleanup() {
  local status=$?
  if [[ "${KEEP_TMP:-0}" == "1" ]]; then
    printf 'KEEP_TMP=1 set; keeping temporary directory: %s\n' "$tmpdir" >&2
  else
    rm -rf -- "$tmpdir" 2>/dev/null || true
  fi
  exit "$status"
}
trap cleanup EXIT INT TERM

mktemp_wrap() {
  local f
  f=$(mktemp "${tmpdir}/5end.XXXXXXXX") || {
    printf 'Error: mktemp failed in %s\n' "$tmpdir" >&2
    exit 1
  }
  printf '%s' "$f"
}

printf 'Temporary directory: %s\n' "$tmpdir" >&2
printf 'Alignment filter: exclude flags 2308 / 0x904; MAPQ >= %s\n' "$mapq" >&2

# Two-column chrom sizes: contig \t length.
# bedtools uses this to align bins during unionbedg.
genome_size_file=$(mktemp_wrap)
cut -f1,2 "${genome}.fai" > "$genome_size_file"

# ---------- Strand selection ----------
do_plus=false
do_minus=false
[[ "$strand_choice" == "both" || "$strand_choice" == "plus"  ]] && do_plus=true
[[ "$strand_choice" == "both" || "$strand_choice" == "minus" ]] && do_minus=true

# ---------- Per-sample containers ----------
plus_files=()     # per-sample bedGraph with 5′ counts on '+' strand
minus_files=()    # per-sample bedGraph with 5′ counts on '-' strand
sample_names=()   # clean headers for columns

# Strict primary mapped only:
#   0x4   unmapped
#   0x100 secondary
#   0x800 supplementary
# Decimal total = 4 + 256 + 2048 = 2308.
FILTER_FLAGS=2308

# ---------- Helper: sanitize sample name for header ----------
# We strip extensions and replace whitespace with underscores so headers
# are TSV-friendly and reproducible. If duplicate names arise, we deduplicate.
sanitize_sample() {
  local s="$1"
  s="${s%.bam}"
  s="${s%.sam}"
  s="${s%.cram}"
  s="${s//$'\t'/_}"
  s="${s// /_}"
  printf '%s' "$s"
}

declare -A seen_names

dedup_name() {
  local n="$1"
  local i=1
  while [[ -n "${seen_names[$n]:-}" ]]; do
    n="${n}_$i"
    i=$((i + 1))
  done
  seen_names["$n"]=1
  printf '%s' "$n"
}

# ---------- Core loop: process each BAM ----------
# For each input BAM:
#   1) coordinate sort + index
#   2) filter to strict primary-only plus MAPQ threshold + index
#   3) compute per-base 5′ counts with bedtools genomecov -dz, strand-specific
#   4) format as 1-bp bedGraph rows: chrom  start0  end0+1  count
for bam in "${bam_files[@]}"; do
  [[ -f "$bam" ]] || {
    printf 'Error: BAM not found: %s\n' "$bam" >&2
    exit 1
  }

  # Build a stable, readable sample name.
  base=$(basename "$bam")
  raw=$(sanitize_sample "$base")
  name=$(dedup_name "$raw")
  sample_names+=("$name")

  printf 'Processing sample: %s\n' "$name" >&2

  # (1) Coordinate sort + index.
  # This is required for fast/random access and consistent downstream behavior.
  sorted_bam=$(mktemp_wrap)
  samtools sort -@ "$threads" \
    -T "${tmpdir}/${name}.sorttmp" \
    -o "$sorted_bam" \
    "$bam"

  samtools index -@ "$threads" "$sorted_bam"

  # (2) Filter alignments + index.
  # Keep:
  #   mapped
  #   primary
  #   non-supplementary
  #   MAPQ >= $mapq
  #
  # Note:
  #   --mapq 0 effectively disables MAPQ filtering.
  filtered_bam=$(mktemp_wrap)
  samtools view -@ "$threads" \
    -b \
    -F "$FILTER_FLAGS" \
    -q "$mapq" \
    "$sorted_bam" \
    -o "$filtered_bam"

  samtools index -@ "$threads" "$filtered_bam"

  # (3) & (4) Compute 5′ counts → 1-bp bedGraph.
  #     • -5 tells genomecov to count only the 5′ position of each read.
  #     • -dz emits a line for every base with nonzero count: chrom pos0 count.
  #     • We convert to bedGraph columns: chrom  start0  end0+1  count.
  if $do_plus; then
    plus_bg=$(mktemp_wrap)
    bedtools genomecov -ibam "$filtered_bam" -strand + -5 -dz \
      | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1, $3}' \
      | LC_ALL=C sort -T "$tmpdir" -k1,1 -k2,2n > "$plus_bg"

    plus_files+=("$plus_bg")
  fi

  if $do_minus; then
    minus_bg=$(mktemp_wrap)
    bedtools genomecov -ibam "$filtered_bam" -strand - -5 -dz \
      | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1, $3}' \
      | LC_ALL=C sort -T "$tmpdir" -k1,1 -k2,2n > "$minus_bg"

    minus_files+=("$minus_bg")
  fi
done

# ---------- Helper: merge a strand across samples ----------
# We need to align per-sample bedGraphs so each row has counts from all samples
# at the same genomic bin. bedtools unionbedg does this, but it requires 2+
# inputs. If only one sample has any signal on a strand, we bypass unionbedg
# and expand that sample's counts into an N-column row, filling zeroes elsewhere.
emit_block() {
  # $1 = '+' or '-'
  # remaining arguments = per-sample bedGraphs for that strand
  local strand="$1"
  shift

  local -a files=("$@")
  local N=${#files[@]}   # number of samples/columns

  # Collect only non-empty files.
  # Some samples may have no signal for this strand.
  local -a nonempty=()
  local -a idx=()
  local i

  for i in "${!files[@]}"; do
    if [[ -s "${files[$i]}" ]]; then
      nonempty+=("${files[$i]}")
      idx+=("$((i + 1))")  # sample column index, 1-based
    fi
  done

  # Strand is globally empty across all samples → skip gracefully.
  if [[ ${#nonempty[@]} -eq 0 ]]; then
    printf 'Note: no 5′ signal for strand %s; skipping\n' "$strand" >&2
    return 0
  fi

  local out
  out=$(mktemp_wrap)

  if [[ ${#nonempty[@]} -ge 2 ]]; then
    # Many inputs: align them with unionbedg.
    # unionbedg merges intervals and emits one row per unique bin.
    local merged
    merged=$(mktemp_wrap)

    bedtools unionbedg \
      -i "${nonempty[@]}" \
      -g "$genome_size_file" > "$merged"

    # Map K non-empty columns back into N sample columns, filling zeroes for
    # samples with no data on this strand.
    local idx_csv
    idx_csv=$(IFS=, ; echo "${idx[*]}")

    awk -v N="$N" -v IDX="$idx_csv" -v STR="$strand" '
      BEGIN {
        OFS = "\t"
        split(IDX, keep, ",")
        K = length(keep)
      }
      {
        # unionbedg columns:
        # chrom  start0  end0  v1 v2 ... vK
        chrom = $1
        s0 = $2
        start1 = s0 + 1

        # Reset all sample values to zero.
        for (i = 1; i <= N; i++) {
          v[i] = 0
        }

        # Place each of the K values into the correct N-wide column.
        for (i = 1; i <= K; i++) {
          v[keep[i]] = $(3 + i)
        }

        # Emit:
        # chrom, strand, start1, start1, followed by N counts.
        printf "%s\t%s\t%d\t%d", chrom, STR, start1, start1
        for (i = 1; i <= N; i++) {
          printf "\t%d", v[i]
        }
        print ""
      }
    ' "$merged" > "$out"

  else
    # Exactly one non-empty sample:
    # copy its values into the right column, filling zeroes elsewhere.
    local target="${idx[0]}"

    awk -v N="$N" -v TARGET="$target" -v STR="$strand" '
      BEGIN {
        OFS = "\t"
      }
      {
        # input is bedGraph-like:
        # chrom  start0  end0  val
        chrom = $1
        s0 = $2
        start1 = s0 + 1
        val = $4

        for (i = 1; i <= N; i++) {
          v[i] = 0
        }
        v[TARGET] = val

        printf "%s\t%s\t%d\t%d", chrom, STR, start1, start1
        for (i = 1; i <= N; i++) {
          printf "\t%d", v[i]
        }
        print ""
      }
    ' "${nonempty[0]}" > "$out"
  fi

  printf '%s\n' "$out"
}

# ---------- Build per-strand blocks then combine ----------
out_blocks=()

if $do_plus; then
  b=$(emit_block '+' "${plus_files[@]}")
  [[ -n "${b:-}" ]] && out_blocks+=("$b")
fi

if $do_minus; then
  b=$(emit_block '-' "${minus_files[@]}")
  [[ -n "${b:-}" ]] && out_blocks+=("$b")
fi

# ---------- Emit header ----------
# Always print header, even if no rows end up passing --min-count.
[[ -n "$output" ]] && exec > "$output"

printf 'chromosome\tstrand\tstart\tend'
for s in "${sample_names[@]}"; do
  printf '\t%s' "$s"
done
printf '\n'

# ---------- Sort rows, apply --min-count filter, print ----------
# We sort by chromosome and start, then keep rows whose total across samples
# is >= --min-count.
if [[ ${#out_blocks[@]} -gt 0 ]]; then
  LC_ALL=C sort -T "$tmpdir" -k1,1V -k3,3n "${out_blocks[@]}" \
    | awk -v MIN="$min_count" '
        BEGIN {
          OFS = "\t"
        }
        {
          tot = 0
          for (i = 5; i <= NF; i++) {
            tot += $i
          }
          if (tot >= MIN) {
            print
          }
        }
      '
fi

# All temporary files are cleaned via trap on exit.