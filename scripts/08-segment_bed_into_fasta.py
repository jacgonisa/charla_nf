#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import sys

# Define file paths
fasta_file = sys.argv[1] # "CANDIDATE_reads_thr50.fa"
bed_file = sys.argv[2] #"CANDIDATE_reads_thr50.kmer.bed"
output_file = sys.argv[3] #"segmented_reads.fa"

# Load sequences from the FASTA into a dictionary: {read_id: sequence}
sequences = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}

# Dictionary to keep track of the segment number for each read
segment_counts = defaultdict(int)

with open(bed_file, 'r') as bed, open(output_file, 'w') as out_f:
    for line in bed:
        # Skip empty lines or comments
        if line.startswith("#") or not line.strip():
            continue
        # Split the line (assuming tab-delimited)
        parts = line.strip().split()
        read_id = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        tag = parts[3]

        # Increment segment counter for this read
        segment_counts[read_id] += 1
        segment_number = segment_counts[read_id]

        # Extract the segment from the sequence.
        # NOTE: BED files are typically 0-based and half-open.
        # If your FASTA uses 1-based indexing, adjust the indices accordingly.
        segment_seq = sequences[read_id][start:end]

        # Create new header following the format: oldreadname_segmentnumber_tag
        new_header = f">{read_id}_{segment_number}_{tag}"
        
        # Write to the output FASTA file
        out_f.write(new_header + "\n")
        # Optionally wrap the sequence at 60 characters per line
        for i in range(0, len(segment_seq), 60):
            out_f.write(str(segment_seq[i:i+60]) + "\n")

