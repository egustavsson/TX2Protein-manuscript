#!/bin/bash

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --directory) directory="$2"; shift ;;
        --output) output_directory="$2"; shift ;;
        --gene) gene_id="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate required parameters
if [[ -z $directory || -z $output_directory || -z $gene_id ]]; then
    echo "Error: Missing required parameter(s)."
    echo "Usage: bash script_name.sh --directory /path/to/directory/ --output /path/to/output/ --gene ENSG00000157087.18"
    exit 1
fi

for i in "$directory"*.gtf; do
    output_file="${output_directory}$(basename "$i")"
    cat "$i" | grep -E "$gene_id" > "$output_file"
done

for i in "$directory"*.tsv; do
    output_file="${output_directory}$(basename "$i")"
    head -n 1 "$i" > "$output_file"
    cat "$i" | grep -E "$gene_id" >> "$output_file"
done
