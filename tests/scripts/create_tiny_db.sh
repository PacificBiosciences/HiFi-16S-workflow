#!/usr/bin/env bash

set -euo pipefail

DBS="${1:-${DBS:-}}"

if [[ -z "${DBS}" ]]; then
  echo "Usage: $0 <database_root>" >&2
  echo "   or: DBS=/path/to/databases $0" >&2
  exit 1
fi

DBS="$(realpath "${DBS}")"

EUK_LONG_IDS="tests/assets/euk_long_reference_ids.txt"
SHARED_DB_SCRIPT="tests/scripts/create_shared_euk_test_db.py"

# =============================================================================
# 16S databases
#
# SILVA, GTDB and Greengenes2 currently support Naive Bayes classification only.
# =============================================================================

for db in silva gtdb gg2; do

  echo
  echo "======================================================================"
  echo "Building tiny 16S database: ${db}"
  echo "======================================================================"

  mkdir -p "tests/tiny_db/${db}/nb"

  nb_src=$(find "${DBS}/${db}/nb" \
    -maxdepth 1 \
    -name '*.fa.gz' \
    -print \
    -quit)

  if [[ -z "${nb_src}" ]]; then
    echo "ERROR: No NB database found for ${db}" >&2
    exit 1
  fi

  nb_name=$(basename "${nb_src}")

  gunzip -c "${nb_src}" |
    awk '/^>/{n++} n<=100' |
    gzip -c \
      >"tests/tiny_db/${db}/nb/${nb_name}"

  gzip -t "tests/tiny_db/${db}/nb/${nb_name}"

  nb_count=$(
    zgrep -c '^>' \
      "tests/tiny_db/${db}/nb/${nb_name}" ||
      true
  )

  echo "NB references: ${nb_count}"

done

# =============================================================================
# EUKARYOME databases
#
# These contain both:
#   - DADA2 / Naive Bayes database
#   - VSEARCH / SINTAX database
# =============================================================================

for db in euk_ssu euk_lsu euk_its euk_long; do

  echo
  echo "======================================================================"
  echo "Building tiny EUKARYOME database: ${db}"
  echo "======================================================================"

  mkdir -p \
    "tests/tiny_db/${db}/nb" \
    "tests/tiny_db/${db}/vsearch"

  # ---------------------------------------------------------------------------
  # Locate source databases
  # ---------------------------------------------------------------------------

  nb_src=$(find "${DBS}/${db}/nb" \
    -maxdepth 1 \
    -name '*.fa.gz' \
    -print \
    -quit)

  vsearch_src=$(find "${DBS}/${db}/vsearch" \
    -maxdepth 1 \
    -name '*.fa.gz' \
    -print \
    -quit)

  if [[ -z "${nb_src}" ]]; then
    echo "ERROR: No NB database found for ${db}" >&2
    exit 1
  fi

  if [[ -z "${vsearch_src}" ]]; then
    echo "ERROR: No VSEARCH database found for ${db}" >&2
    exit 1
  fi

  nb_name=$(basename "${nb_src}")
  vsearch_name=$(basename "${vsearch_src}")

  # ===========================================================================
  # EUK LONG
  #
  # Build NB and VSEARCH from exactly the same selected nucleotide sequences.
  # This makes the test a controlled comparison of the classifiers.
  # ===========================================================================

  if [[ "${db}" == "euk_long" ]]; then

    if [[ ! -f "${EUK_LONG_IDS}" ]]; then
      echo "ERROR: Missing reference ID file:" >&2
      echo "  ${EUK_LONG_IDS}" >&2
      exit 1
    fi

    if [[ ! -f "${SHARED_DB_SCRIPT}" ]]; then
      echo "ERROR: Missing helper script:" >&2
      echo "  ${SHARED_DB_SCRIPT}" >&2
      exit 1
    fi

    python "${SHARED_DB_SCRIPT}" \
      --input "${vsearch_src}" \
      --ids "${EUK_LONG_IDS}" \
      --vsearch-output "tests/tiny_db/${db}/vsearch/${vsearch_name}" \
      --nb-output "tests/tiny_db/${db}/nb/${nb_name}"

  # ===========================================================================
  # EUK ITS
  #
  # NB:
  #   generic first 100 references
  #
  # VSEARCH:
  #   curated references giving known-positive behavior in the test data.
  # ===========================================================================

  elif [[ "${db}" == "euk_its" ]]; then

    # -------------------------------------------------------------------------
    # Naive Bayes
    # -------------------------------------------------------------------------

    gunzip -c "${nb_src}" |
      awk '/^>/{n++} n<=100' |
      gzip -c \
        >"tests/tiny_db/${db}/nb/${nb_name}"

    # -------------------------------------------------------------------------
    # VSEARCH / SINTAX
    # -------------------------------------------------------------------------

    gunzip -c "${vsearch_src}" |
      awk '
        BEGIN {
          # Aspergillus fumigatus
          keep["EUK0551660"] = 1
          keep["EUK0660280"] = 1
          keep["EUK0701466"] = 1
          keep["EUK0989341"] = 1
          keep["EUK1353429"] = 1
          keep["KU761136"]   = 1

          # Cutaneotrichosporon
          keep["EUK0558349"] = 1
          keep["KY103004"]   = 1
          keep["KY103007"]   = 1
          keep["KY103013"]   = 1
          keep["MG857732"]   = 1
          keep["AB018031"]   = 1
          keep["GU941375"]   = 1
          keep["GU941392"]   = 1
          keep["AB305103"]   = 1
        }

        /^>/ {
          header = substr($0, 2)
          split(header, fields, ";")
          id = fields[1]

          print_record = (id in keep)
        }

        print_record
      ' |
      gzip -c \
        >"tests/tiny_db/${db}/vsearch/${vsearch_name}"

  # ===========================================================================
  # EUK SSU / LSU
  #
  # Generic first 100 references from each classifier-specific database.
  # ===========================================================================

  else

    gunzip -c "${nb_src}" |
      awk '/^>/{n++} n<=100' |
      gzip -c \
        >"tests/tiny_db/${db}/nb/${nb_name}"

    gunzip -c "${vsearch_src}" |
      awk '/^>/{n++} n<=100' |
      gzip -c \
        >"tests/tiny_db/${db}/vsearch/${vsearch_name}"

  fi

  # ---------------------------------------------------------------------------
  # Validate
  # ---------------------------------------------------------------------------

  gzip -t "tests/tiny_db/${db}/nb/${nb_name}"
  gzip -t "tests/tiny_db/${db}/vsearch/${vsearch_name}"

  nb_count=$(
    zgrep -c '^>' \
      "tests/tiny_db/${db}/nb/${nb_name}" ||
      true
  )

  vsearch_count=$(
    zgrep -c '^>' \
      "tests/tiny_db/${db}/vsearch/${vsearch_name}" ||
      true
  )

  echo "NB references:      ${nb_count}"
  echo "VSEARCH references: ${vsearch_count}"

  # For euk_long these MUST be the same set size.
  if [[ "${db}" == "euk_long" ]]; then

    if [[ "${nb_count}" -ne "${vsearch_count}" ]]; then
      echo "ERROR: euk_long NB and VSEARCH reference counts differ" >&2
      exit 1
    fi

    echo "euk_long matched reference sets: OK"
  fi

done

echo
echo "======================================================================"
echo "All tiny databases rebuilt successfully."
echo "======================================================================"
