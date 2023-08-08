#!/bin/bash

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <output_dir> <tsv_url>"
    exit 1
fi

# User-defined output directory
output_dir="$1"

# URL of the TSV file
tsv_url="$2"

# Filename for the downloaded TSV file
downloaded_file="$output_dir/panel_data.tsv"

# Filename for the extracted gene set TSV file
gene_set_file="$output_dir/gene_set.tsv"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Download the TSV file
wget "$tsv_url" -O "$downloaded_file"

if [ $? -eq 0 ]; then
    echo "Downloaded TSV file successfully."
else
    echo "Failed to download TSV file."
    exit 1
fi

# Extract the desired column and remove the header
awk -F'\t' '{print $3}' "$downloaded_file" | tail -n +2 > "$gene_set_file"

if [ $? -eq 0 ]; then
    echo "Extracted gene set column and saved to gene_set.tsv."
else
    echo "Failed to extract gene set column."
    exit 1
fi
