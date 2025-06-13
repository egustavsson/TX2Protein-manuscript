#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def load_id_set(id_file):
    with open(id_file, "r", encoding="utf-8") as f:
        ids = set(line.strip() for line in f if line.strip())
    logging.info(f"Loaded {len(ids)} transcript IDs.")
    return ids

def process_gtf_files(gtf_files, id_set, output_path):
    matches = 0
    with open(output_path, "w", encoding="utf-8") as out_f:
        for gtf_file in gtf_files:
            logging.info(f"Searching: {gtf_file}")
            with open(gtf_file, "r", encoding="utf-8") as in_f:
                for line in in_f:
                    if any(tid in line for tid in id_set):
                        out_f.write(line)
                        matches += 1
    logging.info(f"Matched {matches} lines written to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Extract lines matching transcript IDs from GTF files.")
    parser.add_argument("id_file", help="File containing transcript IDs, one per line")
    parser.add_argument("output_file", help="Output file for matching GTF lines")
    parser.add_argument("gtf_files", nargs="+", help="List of GTF files to search")
    args = parser.parse_args()

    id_set = load_id_set(args.id_file)
    process_gtf_files(args.gtf_files, id_set, args.output_file)

if __name__ == "__main__":
    main()

