#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse

def parse_folder_classification(folder, pattern="*.paf"):
    """
    Parse all summary files in a folder to extract read classifications.
    It looks for header lines starting with '#' that have the format:
       #<read_id>|<classification>:<segment_info> [additional details...]
    Instead of concatenating duplicate classifications, they are stored uniquely.
    Returns a dictionary mapping read_id to a semicolon-separated string of unique classifications.
    """
    classification = {}
    filepaths = glob.glob(os.path.join(folder, pattern))
    for filepath in filepaths:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    continue
                # Remove the leading '#' and split on the first pipe
                try:
                    header = line[1:]
                    parts = header.split("|", 1)
                    if len(parts) < 2:
                        continue
                    read_id = parts[0].strip()
                    # Extract classification; take the part before any colon
                    class_part = parts[1].strip()
                    classification_type = class_part.split(":", 1)[0].strip()
                    # Use a set to collect unique classification types
                    if read_id not in classification:
                        classification[read_id] = set()
                    classification[read_id].add(classification_type)
                except Exception as e:
                    print(f"Error parsing line: {line}\n{e}")
    # Convert sets to a semicolon-separated string
    for read_id in classification:
        classification[read_id] = ";".join(sorted(classification[read_id]))
    return classification

def build_summary_table(folders):
    """
    Given a list of summary folders, build a summary table (DataFrame)
    where rows are read IDs and columns are summary folder names.
    The cell values are the unique classifications extracted from that folder.
    """
    summary_data = {}
    for folder in folders:
        col_label = os.path.basename(os.path.normpath(folder))
        print(f"Processing folder: {col_label}")
        folder_class = parse_folder_classification(folder)
        for read_id, cls in folder_class.items():
            if read_id not in summary_data:
                summary_data[read_id] = {}
            summary_data[read_id][col_label] = cls
    df = pd.DataFrame.from_dict(summary_data, orient="index")
    df.index.name = "read_id"
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Generate a summary table comparing classifications from different summary folders."
    )
    parser.add_argument("folders", nargs="+", help="Summary folders (each containing summary PAF files)")
    parser.add_argument("-o", "--outfile", default="summary_table.csv", help="Output CSV filename")
    args = parser.parse_args()

    df = build_summary_table(args.folders)
    df.to_csv(args.outfile)
    print(f"Summary table written to {args.outfile}")

if __name__ == "__main__":
    main()

