#!/usr/bin/env bash
# ================================================================
# 2026-06-04 (v1.2) Jeppe Vinther (with ChatGPT)
#
# BigWig generation for CapZyme-seq / 5′-end display
# ---------------------------------------------------
# Produces strand-specific 5′-end BigWigs for IGV/UCSC display.
#
# Filtering policy
# ----------------
# • Strict primary mapped alignments only:
#     samtools view -F 2308
#
#   This excludes:
#     0x4   unmapped
#     0x100 secondary
#     0x800 supplementary
#
# • MAPQ filter:
#     default: --min-mapq 20
#
#   Change with --min-mapq N.
#   Use --min-mapq 0 to disable MAPQ filtering.
#
# • Duplicates are KEPT by default.
#   If you have UMIs, deduplicate upstream before running this script.
#
# • RPM scaling denominator = reads remaining after the same filtering.
#
# • No min-count filtering is applied to BigWigs.
#   BigWigs are for display/QC and should show all passing 5′ ends.
# ================================================================
#
# Usage:
#   make_tss_bigwig_v1_2.sh [--no-rt] [--keep-all] [--threads N] [--min-mapq N] \
#       <input.{bam|sam|cram|bam.gz|sam.gz}> <reference.fa> <prefix>
#
# Examples:
#   # Default: MAPQ >= 20, exclude unmapped/secondary/supplementary
#   make_tss_bigwig_v1_2.sh --threads 8 sample.bam reference.fa sample
#
#   # Disable MAPQ filtering for comparison
#   make_tss_bigwig_v1_2.sh --min-mapq 0 --threads 8 sample.bam reference.fa sample.no_mapq
#
#   # Plus strand only, useful for plus-only RNA-reference mappings
#   make_tss_bigwig_v1_2.sh --no-rt --threads 8 sample.bam reference.fa sample
# ================================================================

set -euo pipefail
LC_ALL=C

NO_RT=false
KEEP_ALL=false
THREADS=1
MIN_MAPQ=20

# Strict primary mapped only:
#   0x4   unmapped
#   0x100 secondary
#   0x800 supplementary
# Decimal total = 4 + 256 + 2048 = 2308.
FILTER_FLAGS=2308

# -------- robust, order-independent parser --------
positional=()
i=1

while [[ $i -le $# ]]; do
  a="${!i}"

  case "$a" in
    --no-rt)
      NO_RT=true
      ;;

    --keep-all)
      KEEP_ALL=true
      ;;

    --threads)
      ((i++))
      [[ $i -le $# ]] || {
        echo "ERROR: --threads needs a value" >&2
        exit 2
      }
      THREADS="${!i}"
      ;;

    --min-mapq)
      ((i++))
      [[ $i -le $# ]] || {
        echo "ERROR: --min-mapq needs a value" >&2
        exit 2
      }
      MIN_MAPQ="${!i}"
      ;;

    -h|--help)
      sed -n '2,220p' "$0"
      exit 0
      ;;

    --)
      # Explicit separator; ignored here because parser is order-independent.
      ;;

    -*)
      echo "ERROR: unknown option: $a" >&2
      exit 2
      ;;

    *)
      positional+=("$a")
      ;;
  esac

  ((i++))
done

if [[ ${#positional[@]} -ne 3 ]]; then
  echo "Usage: $0 [--no-rt] [--keep-all] [--threads N] [--min-mapq N] <in.bam|sam|cram|gz> <ref.fa> <prefix>" >&2
  printf 'Got %d positional args:\n' "${#positional[@]}" >&2
  for p in "${positional[@]}"; do
    printf '  [%s]\n' "$p" >&2
  done
  exit 2
fi

IN="${positional[0]}"
FA="${positional[1]}"
OUT_PREFIX="$(basename "${positional[2]}")"

# -------- validate option values --------
[[ "$THREADS" =~ ^[0-9]+$ ]] && [[ "$THREADS" -ge 1 ]] || {
  echo "ERROR: --threads must be an integer >= 1" >&2
  exit 2
}

[[ "$MIN_MAPQ" =~ ^[0-9]+$ ]] || {
  echo "ERROR: --min-mapq must be a non-negative integer" >&2
  exit 2
}

# -------- sanity checks --------
for tool in samtools bedtools bedGraphToBigWig; do
  command -v "$tool" >/dev/null || {
    echo "ERROR: $tool not in PATH" >&2
    exit 127
  }
done

[[ -f "$IN" ]] || {
  echo "ERROR: input not found: $IN" >&2
  exit 1
}

[[ -f "$FA" ]] || {
  echo "ERROR: FASTA not found: $FA" >&2
  exit 1
}

# -------- temp cleanup unless --keep-all --------
_tmp_to_clean=()

cleanup() {
  if ! $KEEP_ALL; then
    rm -f "${_tmp_to_clean[@]}" 2>/dev/null || true
  fi
}
trap cleanup EXIT

# -------- helpers --------
make_min_bg() {
  # Create a minimal zero bedGraph so BigWig can be built for empty tracks.
  local sizes="$1"
  local out_bg="$2"
  local chrom
  local len

  if read -r chrom len < "$sizes"; then
    printf "%s\t0\t1\t0\n" "$chrom" > "$out_bg"
  else
    : > "$out_bg"
  fi
}

filter_bg_to_sizes() {
  # Keep only contigs present in CHRSIZES; warn once per dropped contig.
  local sizes="$1"

  awk -v sz="$sizes" '
    BEGIN {
      while ((getline < sz) > 0) {
        ok[$1] = 1
      }
    }
    {
      if (ok[$1]) {
        print
      } else if (!warn[$1]++) {
        printf "Note: dropping contig not in sizes: %s\n", $1 > "/dev/stderr"
      }
    }
  '
}

needs_ref_for_view() {
  case "$IN" in
    *.cram) return 0 ;;
    *)      return 1 ;;
  esac
}

# -------- 0) normalize input to BAM --------
RAW_BAM=""

case "$IN" in
  *.bam)
    RAW_BAM="$IN"
    ;;

  *)
    RAW_BAM="$(mktemp "${OUT_PREFIX}.raw.XXXX.bam")"
    _tmp_to_clean+=("$RAW_BAM")

    echo "► Converting $(basename "$IN") → BAM …"

    if needs_ref_for_view; then
      samtools view -@ "$THREADS" -T "$FA" -b "$IN" -o "$RAW_BAM"
    else
      samtools view -@ "$THREADS" -b "$IN" -o "$RAW_BAM"
    fi
    ;;
esac

# -------- 1) sort + index --------
SORTED_BAM="${OUT_PREFIX}.sorted.bam"

echo "► Sorting BAM (threads=$THREADS) …"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$RAW_BAM"
samtools index -@ "$THREADS" "$SORTED_BAM"

# -------- 2) filter alignments --------
FILTERED_BAM="${OUT_PREFIX}.filtered.bam"

echo "► Filtering alignments:"
echo "   exclude flags: $FILTER_FLAGS / 0x904"
echo "   MAPQ >= $MIN_MAPQ"
echo "   duplicates kept"

samtools view -@ "$THREADS" \
  -b \
  -F "$FILTER_FLAGS" \
  -q "$MIN_MAPQ" \
  "$SORTED_BAM" \
  -o "$FILTERED_BAM"

samtools index -@ "$THREADS" "$FILTERED_BAM"

# -------- 3) chrom sizes (.fai) --------
# Rebuild the .fai if missing OR if FASTA is newer than the index.
if [[ ! -f "${FA}.fai" || "$FA" -nt "${FA}.fai" ]]; then
  echo "► (Re)indexing FASTA with samtools faidx …"
  samtools faidx "$FA"
fi

CHRSIZES="${OUT_PREFIX}.chrom.sizes"
cut -f1,2 "${FA}.fai" > "$CHRSIZES"

[[ -s "$CHRSIZES" ]] || {
  echo "ERROR: empty chrom sizes derived from ${FA}.fai" >&2
  exit 1
}

$KEEP_ALL || _tmp_to_clean+=("$CHRSIZES")

# -------- 4) RPM denominator --------
# Important: denominator is the number of reads left after the same filtering
# used for BigWig signal generation.
TOTAL=$(samtools view -c "$FILTERED_BAM")

if [[ "$TOTAL" -eq 0 ]]; then
  echo "ERROR: No reads left after filtering (-F $FILTER_FLAGS, MAPQ>=$MIN_MAPQ)" >&2
  exit 1
fi

SCALE=$(awk -v n="$TOTAL" 'BEGIN { printf "%.6f", 1000000 / n }')

echo "► $TOTAL reads after filtering → scale ×$SCALE for RPM"

# -------- 5) strand-specific 5′-end coverage --------
STRANDS=( '+' '-' )
$NO_RT && STRANDS=( '+' )

BEDGRAPHS=()

for STR in "${STRANDS[@]}"; do
  label=$([[ "$STR" = "+" ]] && echo plus || echo minus)
  OUTBG="${OUT_PREFIX}.${label}.bedGraph"

  echo "► Computing 5′ coverage for strand $STR ($label) …"

  # bedtools genomecov:
  #   -5       count only the 5′ position of each alignment
  #   -bg      output bedGraph intervals
  #   -strand  strand-specific signal
  #   -scale   RPM scaling using reads after filtering
  #
  # Then:
  #   filter to contigs present in the FASTA .fai
  #   sort for bedGraphToBigWig
  if ! bedtools genomecov \
      -5 \
      -bg \
      -strand "$STR" \
      -scale "$SCALE" \
      -ibam "$FILTERED_BAM" \
      | filter_bg_to_sizes "$CHRSIZES" \
      | LC_ALL=C LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n > "$OUTBG"; then

    echo "ERROR: bedtools genomecov failed for strand $STR" >&2
    exit 1
  fi

  # If empty after filtering, inject a minimal zero entry so BigWig builds.
  if [[ ! -s "$OUTBG" ]]; then
    echo "Note: no 5′ signal on '$label' strand after filtering; emitting empty track." >&2
    make_min_bg "$CHRSIZES" "$OUTBG"
  fi

  # Extra safety: use UCSC bedSort if present.
  if command -v bedSort >/dev/null 2>&1; then
    tmp_sorted="$(mktemp "${OUT_PREFIX}.${label}.XXXX.bg")"
    bedSort "$OUTBG" "$tmp_sorted" && mv -f "$tmp_sorted" "$OUTBG"
  fi

  BEDGRAPHS+=("$OUTBG")

  echo "► Converting $OUTBG → BigWig …"
  bedGraphToBigWig "$OUTBG" "$CHRSIZES" "${OUTBG/.bedGraph/.bw}"

  $KEEP_ALL || _tmp_to_clean+=("$OUTBG")
done

# -------- 6) cleanup / reporting --------
if $KEEP_ALL; then
  echo "► Done. Outputs in: $(pwd)"
  echo "   ${OUT_PREFIX}.plus.bw"
  $NO_RT || echo "   ${OUT_PREFIX}.minus.bw"
  echo "   $SORTED_BAM and ${SORTED_BAM}.bai"
  echo "   $FILTERED_BAM and ${FILTERED_BAM}.bai"
  for bg in "${BEDGRAPHS[@]}"; do
    echo "   $bg"
  done
  echo "   $CHRSIZES"
else
  rm -f "$SORTED_BAM" "${SORTED_BAM}.bai" "$FILTERED_BAM" "${FILTERED_BAM}.bai"

  echo "► Done. Kept only BigWigs in: $(pwd)"
  echo "   ${OUT_PREFIX}.plus.bw"
  $NO_RT || echo "   ${OUT_PREFIX}.minus.bw"
fi