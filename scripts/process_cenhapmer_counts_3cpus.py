#!/usr/bin/env python3

import os
import glob
import argparse
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from tqdm import tqdm

# -----------------------
# Parse command-line args
# -----------------------
parser = argparse.ArgumentParser(description="Process *.cambridge files into k-mer count matrix")
parser.add_argument('--base_dir', type=str, required=True, help='Directory with *.cambridge files')
parser.add_argument('--output', type=str, required=True, help='Path to output TSV file')
parser.add_argument('--total_reads', type=int, required=True, help='Total number of reads in input file')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs allocated by Nextflow for this process') # New argument
args = parser.parse_args()

base_dir = args.base_dir
output_file = args.output
total_reads_per_file = args.total_reads
allocated_cpus = args.cpus # Use the allocated CPUs from Nextflow

output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# -----------------------
# Discover input files
# -----------------------
files = glob.glob(os.path.join(base_dir, "**", "*.txt"), recursive=True) # Changed from *_k41.cambridge to *.txt
if not files:
    raise FileNotFoundError(f"No *.txt files found in {base_dir}") # Updated error message
print(f"[INFO] Found {len(files)} *.txt files in: {base_dir}") # Updated info message

# Sorting logic: Col before Ler, CEN before ARMS, sorted by chromosome number
def custom_sort(files):
    order = {"Col": 1, "Ler": 2, "ARMS": 3, "CEN": 4}
    def sort_key(path):
        fname = os.path.basename(path)
        acc = 'Col' if 'Col' in fname else 'Ler' if 'Ler' in fname else ''
        reg = 'CEN' if 'CEN' in fname else 'ARMS' if 'ARMS' in fname else ''
        chr_num = ''.join(c for c in fname if c.isdigit())
        return (order.get(acc, 0), order.get(reg, 0), chr_num)
    return sorted(files, key=sort_key)

files_sorted = custom_sort(files)
print(f"[INFO] Sorted files: {[os.path.basename(f) for f in files_sorted]}")

# -----------------------
# Processing logic
# -----------------------
reads = defaultdict(dict)

def process_counts_chunk(chunk):
    results = []
    for header, count_line in chunk:
        read_name = header.strip().split()[0][1:]
        count_line = count_line.strip().strip(',')
        counts = np.array(count_line.split(','), dtype=int)
        non_zero = np.count_nonzero(counts)
        results.append((read_name, non_zero))
    return results

# Determine internal parallelism based on allocated_cpus
# If we have multiple files, we'll use ThreadPoolExecutor for files,
# and ProcessPoolExecutor for chunks within each file.
# A simple strategy: use most CPUs for file-level parallelism,
# and 1 CPU for chunk-level parallelism within each file.
num_files_to_process_concurrently = min(allocated_cpus, len(files_sorted))
cpus_per_file_chunk = max(1, allocated_cpus // num_files_to_process_concurrently) # Ensure at least 1 CPU per chunk processing

def process_file(file_path, total_reads, progress_bar, batch_size=5000, cpus_for_chunks=1):
    file_reads = {}
    with open(file_path, 'r') as f:
        batch = []
        header = None
        for line in f:
            if line.startswith(">"):
                header = line
            else:
                if header:
                    batch.append((header, line))
                    header = None
                if len(batch) >= batch_size:
                    with ProcessPoolExecutor(max_workers=cpus_for_chunks) as pool: # Use cpus_for_chunks here
                        result = pool.submit(process_counts_chunk, batch).result()
                        for read, count in result:
                            file_reads[read] = count
                    progress_bar.update(len(batch))
                    batch = []
        if batch:
            with ProcessPoolExecutor(max_workers=cpus_for_chunks) as pool: # Use cpus_for_chunks here
                result = pool.submit(process_counts_chunk, batch).result()
                for read, count in result:
                    file_reads[read] = count
            progress_bar.update(len(batch))
    progress_bar.close()
    return os.path.basename(file_path), file_reads

print(f"[INFO] Using {num_files_to_process_concurrently} threads for parallel file processing, and {cpus_per_file_chunk} CPU(s) per file chunk processing.")

with ThreadPoolExecutor(max_workers=num_files_to_process_concurrently) as executor: # Use num_files_to_process_concurrently
    progress_bars = {
        f: tqdm(
            total=total_reads_per_file,
            desc=os.path.basename(f),
            unit="read",
            leave=False
        ) for f in files_sorted
    }

    futures = {
        executor.submit(process_file, f, total_reads_per_file, progress_bars[f], cpus_for_chunks=cpus_per_file_chunk): f # Pass cpus_per_file_chunk
        for f in files_sorted
    }

    for idx, future in enumerate(as_completed(futures), 1):
        file_path = futures[future]
        try:
            fname, file_reads = future.result()
            for read, count in file_reads.items():
                reads[read][fname] = count
            print(f"[✓] Processed: {fname}")
        except Exception as e:
            print(f"[✗] Error processing {file_path}: {e}")
        finally:
            progress_bars[file_path].close()

# -----------------------
# Write output table
# -----------------------
header = ["ReadName"] + [os.path.basename(f) for f in files_sorted]
rows = []

for read in sorted(reads.keys()):
    row = [read] + [str(reads[read].get(f, 0)) for f in header[1:]]
    rows.append(row)

with open(output_file, 'w') as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

print(f"[✔] Output written to: {output_file}")
