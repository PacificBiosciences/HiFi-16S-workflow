#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.fasta>"
    exit 1
fi

input_file=$1
counter=1
output_fasta="processed.fasta"
output_tsv="output.tsv"

# Clear or create output files
> "$output_fasta"
> "$output_tsv"

# Process the FASTA file
while read line; do
    # Check if line starts with '>' (header line)
    if [[ $line =~ ^\> ]]; then
        # Remove the '>' character and store header
        header=${line#>}
        # Create sequence name
        seqname="seq${counter}"
        # Output to TSV file
        echo -e "${seqname}\t${header}" >> "$output_tsv"
        # Output modified header to FASTA file
        echo ">${seqname} ${header}" >> "$output_fasta"
        ((counter++))
    else
        # Output sequence lines to FASTA file unchanged
        echo "$line" >> "$output_fasta"
    fi
done < "$input_file"

echo "Created $output_tsv and $output_fasta"