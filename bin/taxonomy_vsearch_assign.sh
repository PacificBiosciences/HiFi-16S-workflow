#!/usr/bin/env bash
set -euo pipefail

ASV_FASTA="$1"
DB_FASTA="$2"
DB_TAX="$3"
DB_NAME="$4"
THREADS="$5"
MAXREJECT="$6"
MAXACCEPT="$7"
IDENTITY="$8"
OUT_TSV="$9"

BLAST6="${DB_NAME}.vsearch.blast6"

vsearch \
  --usearch_global "${ASV_FASTA}" \
  --db "${DB_FASTA}" \
  --id "${IDENTITY}" \
  --strand both \
  --threads "${THREADS}" \
  --maxaccepts "${MAXACCEPT}" \
  --maxrejects "${MAXREJECT}" \
  --top_hits_only \
  --blast6out "${BLAST6}"

awk '
BEGIN {
  FS="[ \t]+"
  OFS="\t"

  print "Feature ID", \
        "Reference ID", \
        "Identity", \
        "AlignmentLength", \
        "Taxon"
}

FNR==NR {
  id=$1
  $1=""
  sub(/^[ \t]+/, "", $0)

  # Preserve whitespace within taxon names as spaces.
  gsub(/\t/, " ", $0)

  # Remove optional trailing numeric confidence/weight field.
  sub(/[ \t]+[0-9.]+$/, "", $0)

  tax[id]=$0
  next
}

{
  asv=$1
  ref=$2
  pid=$3
  aln_len=$4

  taxonomy = (ref in tax) ? tax[ref] : "Unclassified"

  print asv, ref, pid, aln_len, taxonomy
}
' "${DB_TAX}" "${BLAST6}" >"${OUT_TSV}"
