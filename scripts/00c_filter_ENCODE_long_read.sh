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
    echo "Usage: bash 00cfilter_ENCODE_long_read.sh --directory /path/to/directory/ --output /path/to/output/ --gene-file genes.tsv"
    exit 1
fi

# Read gene IDs from the gene file into an array
mapfile -t gene_ids < <(grep -v '^$' "$gene_file")  # Remove blank lines

# Generate a pattern for grep using gene IDs
grep_pattern=$(IFS="|"; echo "${gene_ids[*]}")

# Process GTF files
for i in "$directory"*.gtf; do
    output_file="${output_directory}$(basename "$i")"
    grep -E "$grep_pattern" "$i" > "$output_file"
done

# Process TSV files
for i in "$directory"*.tsv; do
    output_file="${output_directory}$(basename "$i")"
    head -n 1 "$i" > "$output_file"
    grep -E "$grep_pattern" "$i" >> "$output_file"
done
