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
args = parser.parse_args()

base_dir = args.base_dir
output_file = args.output
total_reads_per_file = args.total_reads

output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# -----------------------
# Discover input files
# -----------------------
files = glob.glob(os.path.join(base_dir, "**", "*_k41.cambridge"), recursive=True)
if not files:
    raise FileNotFoundError(f"No *_k41.cambridge files found in {base_dir}")
print(f"[INFO] Found {len(files)} *.cambridge files in: {base_dir}")

# Sorting logic: Col before Ler, CEN before ARMS, sorted by chromosome number
def custom_sort(files):
    order = {"Col": 1, "Ler": 2, "ARMS": 3, "CEN": 4}
    def sort_key(path):
        fname = os.path.basename(path)
        acc = 'Col' if 'Col' in fname else 'Ler' if 'Ler' in fname else ''
        reg = 'CEN' if 'CEN' in fname else 'ARMS' if 'ARMS' in fname else ''
        chrn = ''.join(c for c in fname if c.isdigit())
        return (order.get(acc, 0), order.get(reg, 0), chrn)
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

def process_file(file_path, total_reads, progress_bar, batch_size=5000, cpus_per_file=3):
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
                    with ProcessPoolExecutor(max_workers=cpus_per_file) as pool:
                        result = pool.submit(process_counts_chunk, batch).result()
                        for read, count in result:
                            file_reads[read] = count
                    progress_bar.update(len(batch))
                    batch = []
        if batch:
            with ProcessPoolExecutor(max_workers=cpus_per_file) as pool:
                result = pool.submit(process_counts_chunk, batch).result()
                for read, count in result:
                    file_reads[read] = count
            progress_bar.update(len(batch))
    progress_bar.close()
    return os.path.basename(file_path), file_reads

print(f"[INFO] Using 3 CPUs per file, and up to {min(20, len(files_sorted))} threads for parallel processing")
cpus_per_file = 3
num_threads = min(20, len(files_sorted))

with ThreadPoolExecutor(max_workers=num_threads) as executor:
    progress_bars = {
        f: tqdm(
            total=total_reads_per_file,
            desc=os.path.basename(f),
            unit="read",
            leave=False
        ) for f in files_sorted
    }

    futures = {
        executor.submit(process_file, f, total_reads_per_file, progress_bars[f], cpus_per_file=cpus_per_file): f
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

