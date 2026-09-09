#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: prepare_eukaryome_db.sh <db_name> <nb_url> <nb_output> <vsearch_url> <vsearch_output>" >&2
  exit 1
fi

db_name="$1"
nb_url="$2"
nb_output="$3"
vsearch_url="$4"
vsearch_output="$5"

mkdir -p "$(dirname "$nb_output")"
mkdir -p "$(dirname "$vsearch_output")"

# -----------------------------------------------------------------------------
# DADA2 database
# -----------------------------------------------------------------------------

echo "Downloading ${db_name} DADA2 database..."

wget \
  -O dada2.zip \
  "$nb_url"

mkdir dada2_extracted

7z x \
  -y \
  dada2.zip \
  -odada2_extracted

dada2_fasta="$(find dada2_extracted -type f -name '*.fasta' -print -quit)"

# Some EUKARYOME DADA2 downloads contain a nested .7z archive,
# e.g. the long-read database.
if [[ -z "$dada2_fasta" ]]; then

  dada2_7z="$(find dada2_extracted -type f -name '*.7z' -print -quit)"

  if [[ -z "$dada2_7z" ]]; then
    echo "ERROR: no DADA2 FASTA or nested .7z archive found" >&2
    find dada2_extracted -type f -print >&2
    exit 1
  fi

  echo "Found nested DADA2 archive:"
  echo "  $dada2_7z"

  mkdir dada2_inner

  7z x \
    -y \
    "$dada2_7z" \
    -odada2_inner

  dada2_fasta="$(find dada2_inner -type f -name '*.fasta' -print -quit)"
fi

if [[ -z "$dada2_fasta" ]]; then
  echo "ERROR: no DADA2 FASTA found" >&2
  find dada2_extracted -type f -print >&2
  find dada2_inner -type f -print >&2 2>/dev/null || true
  exit 1
fi

echo "Found DADA2 database:"
echo "  $dada2_fasta"

gzip -c "$dada2_fasta" >"$nb_output"

# -----------------------------------------------------------------------------
# SINTAX / VSEARCH database
# -----------------------------------------------------------------------------

echo "Downloading ${db_name} SINTAX database..."

wget \
  -O sintax.zip \
  "$vsearch_url"

mkdir sintax_zip

7z x \
  -y \
  sintax.zip \
  -osintax_zip

sintax_7z="$(find sintax_zip -type f -name '*.7z' -print -quit)"

if [[ -z "$sintax_7z" ]]; then
  echo "ERROR: no .7z archive found inside SINTAX ZIP" >&2
  find sintax_zip -type f -print >&2
  exit 1
fi

echo "Found SINTAX archive:"
echo "  $sintax_7z"

mkdir sintax_extracted

7z x \
  -y \
  "$sintax_7z" \
  -osintax_extracted

sintax_fasta="$(find sintax_extracted -type f -name '*.fasta' -print -quit)"

if [[ -z "$sintax_fasta" ]]; then
  echo "ERROR: no SINTAX FASTA found" >&2
  find sintax_extracted -type f -print >&2
  exit 1
fi

echo "Found SINTAX database:"
echo "  $sintax_fasta"

gzip -c "$sintax_fasta" >"$vsearch_output"

# -----------------------------------------------------------------------------
# Minimal validation
# -----------------------------------------------------------------------------

if ! head -n 1 "$dada2_fasta" | grep -q '^>'; then
  echo "ERROR: DADA2 input is not FASTA" >&2
  exit 1
fi

if ! head -n 1 "$sintax_fasta" | grep -q '^>'; then
  echo "ERROR: SINTAX input is not FASTA" >&2
  exit 1
fi

if ! gzip -t "$nb_output"; then
  echo "ERROR: DADA2 output is not valid gzip" >&2
  exit 1
fi

if ! gzip -t "$vsearch_output"; then
  echo "ERROR: VSEARCH output is not valid gzip" >&2
  exit 1
fi

echo "Finished preparing ${db_name}"
echo "DADA2:   ${nb_output}"
echo "VSEARCH: ${vsearch_output}"
