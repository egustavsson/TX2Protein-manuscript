#!/usr/bin/env python3

import csv
import argparse
import logging
import os
import sys

csv.field_size_limit(sys.maxsize)

MAX_EXCEL_CELL_LENGTH = 16000
EXCEL_ESCAPE_PREFIX = "'"

# Canonical order of shared columns
SHARED_COLUMNS = [
    "Protein Accession", "Gene", "Organism", "Protein Full Name", "Protein Unmodified Mass",
    "Number of Proteins in Group", "Unique Peptides", "Shared Peptides", "Number of Peptides",
    "Number of Unique Peptides", "Sequence Coverage Fraction", "Sequence Coverage",
    "Sequence Coverage with Mods", "Modification Info List", "Number of PSMs",
    "Protein Decoy/Contaminant/Target", "Protein Cumulative Target",
    "Protein Cumulative Decoy", "Protein QValue", "Best Peptide Score",
    "Best Peptide Notch QValue"
]

def load_valid_ids(path):
    with open(path, "r", encoding="utf-8") as f:
        return set(line.strip() for line in f if line.strip())

def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    return header, rows

def extract_column_indices(header):
    return {
        'protein': header.index("Protein Accession"),
        'shared': header.index("Shared Peptides"),
        'unique': header.index("Unique Peptides"),
        'qval': header.index("Protein QValue"),
    }

def partition_rows(rows, protein_idx, valid_ids):
    included, excluded = [], []
    for row in rows:
        proteins = row[protein_idx].strip().split("|")
        if all(pid in valid_ids for pid in proteins):
            included.append(row)
        else:
            excluded.append(row)
    return included, excluded

def collect_excluded_shared_peptides(excluded, shared_idx):
    peptides = set()
    for row in excluded:
        if len(row) > shared_idx:
            for pep in row[shared_idx].strip().split("|"):
                if pep:
                    peptides.add(pep)
    return peptides

def escape_excel(value):
    if not value:
        return "N/A"
    if isinstance(value, str) and value[0] in ('=', '+', '-', '@'):
        return EXCEL_ESCAPE_PREFIX + value
    return value

def truncate_and_escape_fields(row):
    return [
        escape_excel(field[:MAX_EXCEL_CELL_LENGTH] if isinstance(field, str) and len(field) > MAX_EXCEL_CELL_LENGTH else field)
        if field.strip() else "N/A"
        for field in row
    ]

def build_merged_header(all_headers):
    all_columns = set()
    for hdr in all_headers:
        all_columns.update(h.strip() for h in hdr)
    non_shared = sorted(all_columns - set(SHARED_COLUMNS))
    return ["Origin"] + [col for col in SHARED_COLUMNS if col in all_columns] + non_shared

def normalize_row(row_dict, full_header):
    return [row_dict.get(col, "N/A") if row_dict.get(col, "").strip() else "N/A" for col in full_header]

def clean_and_transform_rows(included, header, shared_idx, unique_idx, qval_idx, excluded_peptides, protein_idx, origin_path):
    header_map = {h.strip(): i for i, h in enumerate(header)}
    cleaned = []
    for row in included:
        shared = row[shared_idx].strip().split("|") if row[shared_idx].strip() else []
        unique = row[unique_idx].strip().split("|") if row[unique_idx].strip() else []
        shared = [p for p in shared if p not in excluded_peptides]
        if not unique and not shared:
            continue
        row[shared_idx] = "|".join(shared)
        row[unique_idx] = "|".join(unique)
        try:
            if float(row[qval_idx]) > 0.01:
                continue
        except ValueError:
            continue
        cleaned_row = truncate_and_escape_fields(row)
        row_dict = {col: cleaned_row[idx] for col, idx in header_map.items()}
        row_dict["Origin"] = escape_excel(origin_path)
        cleaned.append(row_dict)
    return cleaned

def write_output(output_path, full_header, all_rows):
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(full_header)
        for row_dict in all_rows:
            writer.writerow(normalize_row(row_dict, full_header))

def main(id_file, output_file, tsv_files):
    logging.basicConfig(level=logging.INFO)
    valid_ids = load_valid_ids(id_file)
    all_rows = []
    all_headers = []

    for tsv_path in tsv_files:
        logging.info(f"Processing: {tsv_path}")
        header, rows = read_tsv(tsv_path)
        all_headers.append(header)
        indices = extract_column_indices(header)
        included, excluded = partition_rows(rows, indices['protein'], valid_ids)
        excluded_peptides = collect_excluded_shared_peptides(excluded, indices['shared'])
        cleaned_dict_rows = clean_and_transform_rows(
            included,
            header,
            indices['shared'],
            indices['unique'],
            indices['qval'],
            excluded_peptides,
            indices['protein'],
            os.path.abspath(tsv_path)
        )
        all_rows.extend(cleaned_dict_rows)

    if all_rows:
        merged_header = build_merged_header(all_headers)
        write_output(output_file, merged_header, all_rows)
        logging.info(f"Final output written to: {output_file}")
    else:
        logging.warning("No valid data to write.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter protein TSV files with ID list.")
    parser.add_argument("id_file", help="Path to unique_transcript_IDs_novel.txt")
    parser.add_argument("output_file", help="Combined output TSV file path")
    parser.add_argument("tsv_files", nargs="+", help="One or more AllQuantifiedProteinGroups.tsv files")
    args = parser.parse_args()
    main(args.id_file, args.output_file, args.tsv_files)
