#!/bin/bash

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --directory) directory="$2"; shift ;;
        --output) output_directory="$2"; shift ;;
        --gene-file) gene_file="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate required parameters
if [[ -z $directory || -z $output_directory || -z $gene_file ]]; then
    echo "Error: Missing required parameter(s)."
    echo "Usage: bash script_name.sh --directory /path/to/directory/ --output /path/to/output/ --gene-file genes.tsv"
    exit 1
fi

# Create a temporary file for storing the gene IDs
temp_gene_file=$(mktemp)

# Extract gene IDs from the TSV file and save them to the temporary file
awk -F'\t' '{print $1}' "$gene_file" > "$temp_gene_file"

for i in "$directory"*.gtf; do
    output_file="${output_directory}$(basename "$i")"
    grep -Ff "$temp_gene_file" "$i" > "$output_file"
done

for i in "$directory"*.tsv; do
    output_file="${output_directory}$(basename "$i")"
    head -n 1 "$i" > "$output_file"
    grep -Ff "$temp_gene_file" "$i" >> "$output_file"
done

# Clean up the temporary gene file
rm "$temp_gene_file"
