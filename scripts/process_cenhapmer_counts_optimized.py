#!/usr/bin/env python3

import os
import glob
import argparse
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, TimeoutError
from tqdm import tqdm
import sys
from multiprocessing import cpu_count
import math

# -----------------------
# Parse command-line args
# -----------------------
parser = argparse.ArgumentParser(description="Process *.txt files into k-mer count matrix")
parser.add_argument('--base_dir', type=str, required=True, help='Directory with *.txt files')
parser.add_argument('--output', type=str, required=True, help='Path to output TSV file')
parser.add_argument('--total_reads', type=int, required=True, help='Total number of reads in input file')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs allocated by Nextflow for this process')
args = parser.parse_args()

base_dir = args.base_dir
output_file = args.output
total_reads_per_file = args.total_reads
allocated_cpus = min(args.cpus, cpu_count())  # Don't exceed system CPU count

output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# -----------------------
# Discover input files
# -----------------------
files = glob.glob(os.path.join(base_dir, "**", "*.txt"), recursive=True)
if not files:
    raise FileNotFoundError(f"No *.txt files found in {base_dir}")
print(f"[INFO] Found {len(files)} *.txt files in: {base_dir}")

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
# Optimized processing logic
# -----------------------
reads = defaultdict(dict)

def process_counts_chunk(chunk):
    """Process a chunk of header/count pairs efficiently"""
    results = []
    for header, count_line in chunk:
        try:
            read_name = header.strip().split()[0][1:]
            count_line = count_line.strip().strip(',')
            counts = np.array(count_line.split(','), dtype=int)
            non_zero = np.count_nonzero(counts)
            results.append((read_name, non_zero))
        except Exception as e:
            print(f"[WARNING] Skipping malformed entry: {e}", file=sys.stderr)
            continue
    return results

# Optimized parallelization strategy
num_files = len(files_sorted)
if allocated_cpus <= 2 or num_files == 1:
    # Use all CPUs for single file processing
    file_workers = 1
    chunk_workers = allocated_cpus
elif num_files >= allocated_cpus:
    # Use file-level parallelism primarily
    file_workers = allocated_cpus
    chunk_workers = 1
else:
    # Balance between file and chunk parallelism
    file_workers = min(num_files, max(1, allocated_cpus // 2))
    chunk_workers = max(1, allocated_cpus // file_workers)

print(f"[INFO] Using {file_workers} file workers and {chunk_workers} chunk workers per file")

# Create a single shared ProcessPoolExecutor for all chunk processing
chunk_executor = ProcessPoolExecutor(max_workers=chunk_workers)

def process_file_optimized(file_path, total_reads, batch_size=20000):
    """Process a single file with optimized batching and shared process pool"""
    print(f"[INFO] Starting to process {os.path.basename(file_path)}")
    
    file_reads = {}
    processed_reads = 0
    
    # Create progress bar
    pbar = tqdm(
        total=total_reads,
        desc=os.path.basename(file_path),
        unit="read",
        leave=False
    )
    
    try:
        with open(file_path, 'r') as f:
            batch = []
            header = None
            pending_futures = []
            
            for line in f:
                if line.startswith(">"):
                    header = line
                else:
                    if header:
                        batch.append((header, line))
                        header = None
                        
                    if len(batch) >= batch_size:
                        # Submit batch for processing
                        future = chunk_executor.submit(process_counts_chunk, batch.copy())
                        pending_futures.append((future, len(batch)))
                        batch = []
                        
                        # Process completed futures to avoid memory buildup
                        completed_futures = []
                        for future, batch_len in pending_futures:
                            if future.done():
                                try:
                                    result = future.result(timeout=0.1)
                                    for read, count in result:
                                        file_reads[read] = count
                                    pbar.update(batch_len)
                                    completed_futures.append((future, batch_len))
                                except TimeoutError:
                                    pass  # Future not ready yet
                                except Exception as e:
                                    print(f"[WARNING] Batch processing failed: {e}", file=sys.stderr)
                                    completed_futures.append((future, batch_len))
                        
                        # Remove completed futures
                        pending_futures = [(f, l) for f, l in pending_futures if (f, l) not in completed_futures]
            
            # Process remaining batch
            if batch:
                future = chunk_executor.submit(process_counts_chunk, batch)
                pending_futures.append((future, len(batch)))
            
            # Wait for all remaining futures
            for future, batch_len in pending_futures:
                try:
                    result = future.result(timeout=300)  # 5 minute timeout per batch
                    for read, count in result:
                        file_reads[read] = count
                    pbar.update(batch_len)
                except TimeoutError:
                    print(f"[WARNING] Batch timed out for {os.path.basename(file_path)}", file=sys.stderr)
                except Exception as e:
                    print(f"[WARNING] Batch processing failed for {os.path.basename(file_path)}: {e}", file=sys.stderr)
        
        pbar.close()
        return os.path.basename(file_path), file_reads
        
    except Exception as e:
        pbar.close()
        print(f"[ERROR] Failed to process {file_path}: {e}", file=sys.stderr)
        return os.path.basename(file_path), {}

# Process files with optimized parallelization
try:
    with ThreadPoolExecutor(max_workers=file_workers) as file_executor:
        # Submit all files for processing
        futures = {
            file_executor.submit(process_file_optimized, f, total_reads_per_file): f 
            for f in files_sorted
        }
        
        # Process completed futures
        for future in as_completed(futures):
            file_path = futures[future]
            try:
                fname, file_reads_data = future.result(timeout=3600)  # 1 hour timeout per file
                
                # Store results
                for read, count in file_reads_data.items():
                    reads[read][fname] = count
                    
                print(f"[✓] Processed: {fname}")
                
            except TimeoutError:
                print(f"[✗] Processing of {os.path.basename(file_path)} timed out after 1 hour", file=sys.stderr)
            except Exception as e:
                print(f"[✗] Error processing {os.path.basename(file_path)}: {e}", file=sys.stderr)
                
except KeyboardInterrupt:
    print("[INFO] Processing interrupted by user")
    sys.exit(1)
finally:
    # Clean up the shared process pool
    chunk_executor.shutdown(wait=True)

# -----------------------
# Write output table
# -----------------------
print("[INFO] Writing output table...")

header = ["ReadName"] + [os.path.basename(f) for f in files_sorted]
rows = []

for read in sorted(reads.keys()):
    row = [read] + [str(reads[read].get(os.path.basename(f), 0)) for f in files_sorted]
    rows.append(row)

try:
    with open(output_file, 'w') as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    
    print(f"[✔] Output written to: {output_file}")
    print(f"[INFO] Processed {len(reads)} unique reads across {len(files_sorted)} files")
    
except Exception as e:
    print(f"[ERROR] Failed to write output: {e}", file=sys.stderr)
    sys.exit(1)
