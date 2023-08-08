#!/bin/bash

# Check if directory argument is provided
if [ $# -eq 0 ]; then
  echo "Usage: ./download_files.sh <directory>"
  exit 1
fi

directory="$1"

# Create the directory if it doesn't exist
mkdir -p "$directory"

# Array of file URLs
urls=(
  "https://www.encodeproject.org/files/ENCFF771XYK/@@download/ENCFF771XYK.tsv"
  "https://www.encodeproject.org/files/ENCFF385SKN/@@download/ENCFF385SKN.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF129ZDE/@@download/ENCFF129ZDE.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF124RDD/@@download/ENCFF124RDD.tsv"
  "https://www.encodeproject.org/files/ENCFF635HXQ/@@download/ENCFF635HXQ.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF107YAW/@@download/ENCFF107YAW.tsv"
  "https://www.encodeproject.org/files/ENCFF509ZIC/@@download/ENCFF509ZIC.tsv"
  "https://www.encodeproject.org/files/ENCFF270HYK/@@download/ENCFF270HYK.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF197DZP/@@download/ENCFF197DZP.tsv"
  "https://www.encodeproject.org/files/ENCFF573QGI/@@download/ENCFF573QGI.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF457BBP/@@download/ENCFF457BBP.tsv"
  "https://www.encodeproject.org/files/ENCFF180NXA/@@download/ENCFF180NXA.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF015EMC/@@download/ENCFF015EMC.tsv"
  "https://www.encodeproject.org/files/ENCFF561RPD/@@download/ENCFF561RPD.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF707NQJ/@@download/ENCFF707NQJ.gtf.gz"
  "https://www.encodeproject.org/files/ENCFF367FEB/@@download/ENCFF367FEB.tsv"
  "https://www.encodeproject.org/files/ENCFF482LTF/@@download/ENCFF482LTF.tsv"
  "https://www.encodeproject.org/files/ENCFF519SVS/@@download/ENCFF519SVS.gtf.gz"
)

total_size=0

# Calculate total file size
for url in "${urls[@]}"; do
  filename=$(basename "$url")
  file_size=$(curl -sI "$url" | awk '/Content-Length/{size=$2} END{print size}')
  file_size="${file_size//[$'\t\r\n']}"
  if [[ $file_size =~ ^[0-9]+$ ]]; then
    total_size=$((total_size + file_size))
  else
    echo "Unable to retrieve file size for $filename"
  fi
done

# Convert total_size to human-readable format
total_size_human=$(numfmt --to=iec --suffix=B "$total_size")

echo "Total file size: $total_size_human"

read -p "Do you want to continue with the download? (y/n): " choice

if [[ "$choice" =~ ^[Yy]$ ]]; then
  # Download files
  for url in "${urls[@]}"; do
    filename=$(basename "$url")
    filepath="$directory/$filename"
    echo "Downloading $filename..."
    curl -# -L "$url" -o "$filepath"
    echo "Downloaded $filename"
  done
else
  echo "Download canceled"
fi
