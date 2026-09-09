#!/usr/bin/env python3

import argparse
import gzip
import re
from pathlib import Path


RANK_ORDER = ["d", "p", "c", "o", "f", "g", "s"]


def open_text(path, mode="rt"):
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, mode)

    return open(path, mode)


def read_fasta(path):
    header = None
    seq = []

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)

                header = line[1:]
                seq = []
            else:
                seq.append(line)

    if header is not None:
        yield header, "".join(seq)


def parse_sintax_header(header):
    """
    Parse a SINTAX-style EUKARYOME header.

    Example:

    EUK1223324;tax=d:Rhizaria,p:Cercozoa,c:Sarcomonadea,
    o:Paracercomonadida,f:Paracercomonadidae,g:Paracercomonas;
    """

    seq_id = header.split(";", 1)[0]

    match = re.search(r";tax=(.*?);?$", header)

    if match is None:
        raise ValueError(f"No SINTAX taxonomy found in header:\n{header}")

    taxonomy_string = match.group(1)

    taxonomy = {}

    for field in taxonomy_string.split(","):
        field = field.strip()

        if not field:
            continue

        if ":" not in field:
            raise ValueError(f"Malformed taxonomy field '{field}' in header:\n{header}")

        rank, value = field.split(":", 1)

        rank = rank.strip()
        value = value.strip()

        if not rank:
            raise ValueError(f"Empty taxonomy rank in header:\n{header}")

        if not value:
            continue

        taxonomy[rank] = value

    return seq_id, taxonomy


def dada2_header(taxonomy):
    """
    Convert SINTAX taxonomy:

        d:Rhizaria,p:Cercozoa,c:Sarcomonadea,...

    into DADA2 taxonomy header format:

        Rhizaria;Cercozoa;Sarcomonadea;...;

    Stop at the first missing rank so that the hierarchy remains contiguous.
    """

    ranks = []

    for rank in RANK_ORDER:
        value = taxonomy.get(rank)

        if value is None or value == "":
            break

        ranks.append(value)

    if not ranks:
        raise ValueError("No usable taxonomy ranks found")

    return ";".join(ranks) + ";"


def read_ids(path):
    ids = set()

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            ids.add(line)

    return ids


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create matched tiny DADA2 and VSEARCH databases "
            "from the same EUKARYOME SINTAX reference sequences."
        )
    )

    parser.add_argument(
        "--input", required=True, help="Full EUKARYOME SINTAX FASTA (.fa or .fa.gz)"
    )

    parser.add_argument(
        "--ids", required=True, help="File containing one reference ID per line"
    )

    parser.add_argument(
        "--vsearch-output",
        required=True,
        help="Output SINTAX FASTA for VSEARCH (.fa.gz)",
    )

    parser.add_argument(
        "--nb-output", required=True, help="Output DADA2 taxonomy FASTA (.fa.gz)"
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    ids_file = Path(args.ids)
    vsearch_output = Path(args.vsearch_output)
    nb_output = Path(args.nb_output)

    if not input_file.exists():
        raise FileNotFoundError(f"Input database not found: {input_file}")

    if not ids_file.exists():
        raise FileNotFoundError(f"ID file not found: {ids_file}")

    wanted = read_ids(ids_file)

    if not wanted:
        raise RuntimeError(f"No sequence IDs found in: {ids_file}")

    vsearch_output.parent.mkdir(parents=True, exist_ok=True)

    nb_output.parent.mkdir(parents=True, exist_ok=True)

    found = set()

    with (
        gzip.open(vsearch_output, "wt") as vsearch_out,
        gzip.open(nb_output, "wt") as nb_out,
    ):
        for header, sequence in read_fasta(input_file):
            # --------------------------------------------------------------
            # Extract only the sequence ID first.
            #
            # This avoids parsing taxonomy for every record in the full
            # EUKARYOME database. Only selected references are parsed.
            # --------------------------------------------------------------

            seq_id = header.split(";", 1)[0]

            if seq_id not in wanted:
                continue

            # --------------------------------------------------------------
            # Parse taxonomy only for selected records.
            # --------------------------------------------------------------

            seq_id, taxonomy = parse_sintax_header(header)

            found.add(seq_id)

            # --------------------------------------------------------------
            # VSEARCH database
            #
            # Keep the original EUKARYOME SINTAX header and exact sequence.
            # --------------------------------------------------------------

            vsearch_out.write(f">{header}\n{sequence}\n")

            # --------------------------------------------------------------
            # DADA2 / Naive Bayes database
            #
            # Same exact nucleotide sequence, but taxonomy is converted to
            # the semicolon-delimited format expected by DADA2.
            # --------------------------------------------------------------

            nb_taxonomy = dada2_header(taxonomy)

            nb_out.write(f">{nb_taxonomy}\n{sequence}\n")

    # ----------------------------------------------------------------------
    # Validate that every requested reference was found.
    # ----------------------------------------------------------------------

    missing = wanted - found

    if missing:
        raise RuntimeError(
            "The following requested IDs were not found "
            "in the SINTAX database:\n" + "\n".join(sorted(missing))
        )

    print(f"Requested references: {len(wanted)}")
    print(f"Selected references:  {len(found)}")
    print()
    print(f"VSEARCH database: {vsearch_output}")
    print(f"DADA2 database:   {nb_output}")


if __name__ == "__main__":
    main()
