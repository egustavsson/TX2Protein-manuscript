#!/usr/bin/env python3

import pandas as pd
import re
import csv
import argparse
from collections import defaultdict

def extract_gtf_annotations(gtf_path):
    transcript_info = {}
    with open(gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or 'transcript_id "' not in line:
                continue
            if "\ttranscript\t" not in line:
                continue  # Only use 'transcript' entries for location

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            attributes = fields[8]

            match = re.search(r'transcript_id "([^"]+)"', attributes)
            if not match:
                continue

            transcript_id = match.group(1)

            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            status_match = re.search(r'transcript_status "([^"]+)"', attributes)

            gene_name = gene_name_match.group(1) if gene_name_match else "N/A"
            transcript_status = status_match.group(1) if status_match else "N/A"

            transcript_info[transcript_id] = {
                "gene_name": gene_name,
                "transcript_status": transcript_status,
                "chromosome": chrom,
                "start": start,
                "end": end
            }
    return transcript_info

def load_cleaned_proteins(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df.fillna("N/A", inplace=True)
    return df

def build_transcript_index(df):
    transcript_to_origins = defaultdict(set)
    transcript_to_peptides = defaultdict(set)

    for _, row in df.iterrows():
        origin = row.get("Origin", "N/A")
        accessions = row.get("Protein Accession", "")
        unique_peps = row.get("Unique Peptides", "")
        shared_peps = row.get("Shared Peptides", "")

        peptides = set()
        if unique_peps and unique_peps != "N/A":
            peptides.update(unique_peps.split("|"))
        if shared_peps and shared_peps != "N/A":
            peptides.update(shared_peps.split("|"))

        for transcript in accessions.split("|"):
            transcript_to_origins[transcript].add(origin)
            transcript_to_peptides[transcript].update(peptides)

    return transcript_to_origins, transcript_to_peptides

def generate_summary(gtf_path, cleaned_tsv_path, output_path):
    gtf_info = extract_gtf_annotations(gtf_path)
    df = load_cleaned_proteins(cleaned_tsv_path)
    transcript_to_origins, transcript_to_peptides = build_transcript_index(df)

    summary = []
    for transcript, origins in transcript_to_origins.items():
        info = gtf_info.get(transcript, {
            "gene_name": "N/A",
            "transcript_status": "N/A",
            "chromosome": "N/A",
            "start": "N/A",
            "end": "N/A"
        })
        peptides = sorted(transcript_to_peptides.get(transcript, set()))
        summary.append([
            transcript,
            info["gene_name"],
            info["transcript_status"],
            len(origins),
            "|".join(peptides) if peptides else "N/A",
            info["chromosome"],
            info["start"],
            info["end"]
        ])

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "transcript_id",
            "gene_name",
            "transcript_status",
            "origin_mass_spec_dataset_count",
            "peptides_all",
            "chromosome",
            "start",
            "end"
        ])
        writer.writerows(summary)

def main():
    parser = argparse.ArgumentParser(description="Summarize transcript info with mass spec support and genomic location.")
    parser.add_argument("gtf_file", help="Path to unique_transcript_IDs_novel_grepped.gtf")
    parser.add_argument("cleaned_tsv", help="Path to cleaned_output.tsv")
    parser.add_argument("output_tsv", help="Path to write summary TSV")
    args = parser.parse_args()

    generate_summary(args.gtf_file, args.cleaned_tsv, args.output_tsv)

if __name__ == "__main__":
    main()
