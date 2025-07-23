#!/usr/bin/env python3
import os
import glob
import argparse
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError # Import TimeoutError
from tqdm import tqdm
# import signal # REMOVE THIS
import sys # Keep sys for sys.exit
# import time # Remove time if not used for artificial pauses

# -----------------------
# Parse command-line args
# -----------------------
parser = argparse.ArgumentParser(description="Process *.txt files into k-mer count matrix")
parser.add_argument('--base_dir', type=str, required=True, help='Directory with *.txt files')
parser.add_argument('--output', type=str, required=True, help='Path to output TSV file')
parser.add_argument('--total_reads', type=int, required=True, help='Total number of reads in input file')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs allocated by Nextflow for this process')
# parser.add_argument('--timeout', type=int, default=3600, help='Timeout per file in seconds (default: 1 hour)') # REMOVE THIS
args = parser.parse_args()

base_dir = args.base_dir
output_file = args.output
total_reads_per_file = args.total_reads
allocated_cpus = args.cpus  # Cap at 8 to avoid over-parallelization (still good practice)
# timeout_per_file = args.timeout # REMOVE THIS

output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# -----------------------
# Timeout handler (REMOVE THIS BLOCK)
# -----------------------
# def timeout_handler(signum, frame):
#     print(f"[ERROR] Process timed out after {timeout_per_file} seconds")
#     sys.exit(1)

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
# Processing logic (simplified, no nested parallelism)
# -----------------------
reads = defaultdict(dict)

def process_counts_chunk(chunk):
    """Process a chunk of header/count pairs - now runs in main thread"""
    results = []
    for header, count_line in chunk:
        try:
            read_name = header.strip().split()[0][1:]
            count_line = count_line.strip().strip(',')
            counts = np.array(count_line.split(','), dtype=int)
            non_zero = np.count_nonzero(counts)
            results.append((read_name, non_zero))
        except Exception as e:
            # Print to stderr for better logging
            print(f"[WARNING] Skipping malformed entry: {e}", file=sys.stderr)
            continue
    return results

def process_file_logic(file_path, total_reads_val, batch_size=10000): # Renamed to avoid confusion with timeout
    """Process a single file without timeout logic here"""
    print(f"[INFO] Starting to process {os.path.basename(file_path)}")
    
    file_reads = {}
    processed_reads = 0
    
    # Create progress bar
    pbar = tqdm(
        total=total_reads_val, # Use total_reads_val passed to the function
        desc=os.path.basename(file_path),
        unit="read",
        leave=False
    )
    
    try:
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
                        result = process_counts_chunk(batch)
                        for read, count in result:
                            file_reads[read] = count
                        
                        processed_reads += len(batch)
                        pbar.update(len(batch))
                        batch = []
                        
                        # No longer need time.sleep(0.001) as it's not yielding for GIL
            
            # Process remaining batch
            if batch:
                result = process_counts_chunk(batch)
                for read, count in result:
                    file_reads[read] = count
                pbar.update(len(batch))
            
        pbar.close()
        # signal.alarm(0) # REMOVE THIS
        
        return os.path.basename(file_path), file_reads
        
    except Exception as e:
        pbar.close() # Ensure pbar is closed on error
        print(f"[ERROR] Failed to process {file_path}: {e}", file=sys.stderr)
        # It's generally better for a threaded function to return an error
        # rather than exiting the whole process directly.
        # The calling ThreadPoolExecutor will catch it.
        # sys.exit(1) # REMOVE THIS
        return os.path.basename(file_path), {} # Return empty on failure for this file

# Use ThreadPoolExecutor for file-level parallelism
# Use 'allocated_cpus' as the max_workers
max_concurrent_files = allocated_cpus # Use allocated_cpus directly
print(f"[INFO] Using {max_concurrent_files} threads for parallel file processing")

# A common timeout for all files, if desired. This will be the maximum time
# any single file processing is allowed to take before being cancelled.
# You could add this as a new command line argument if you want it configurable.
DEFAULT_FILE_TIMEOUT = 3600 # 1 hour
file_processing_timeout = DEFAULT_FILE_TIMEOUT # Can be made a param if needed

try:
    with ThreadPoolExecutor(max_workers=max_concurrent_files) as executor:
        # Submit all files for processing
        futures = {
            executor.submit(process_file_logic, f, total_reads_per_file): f 
            for f in files_sorted
        }
        
        # Process completed futures
        for future in as_completed(futures):
            file_path = futures[future]
            try:
                # Use the timeout parameter here
                fname, file_reads_data = future.result(timeout=file_processing_timeout)
                
                # Store results
                for read, count in file_reads_data.items():
                    reads[read][fname] = count
                    
                print(f"[✓] Processed: {fname}")
                
            except TimeoutError:
                print(f"[✗] Processing of {os.path.basename(file_path)} timed out after {file_processing_timeout} seconds.", file=sys.stderr)
            except Exception as e:
                print(f"[✗] Error processing {os.path.basename(file_path)}: {e}", file=sys.stderr)
                
except KeyboardInterrupt:
    print("[INFO] Processing interrupted by user")
    sys.exit(1)
except Exception as e:
    print(f"[ERROR] Unexpected error during thread pool execution: {e}", file=sys.stderr)
    sys.exit(1)

# -----------------------
# Write output table
# -----------------------
print("[INFO] Writing output table...")

# Ensure header reflects only successfully processed files, or all files if desired
# For robustness, let's build the header from the files_sorted and use get()
# to put 0 for files that failed to process.
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
    print(f"[INFO] Processed {len(reads)} unique reads across {len(files_sorted)} files (considering all found, even if some failed individually)")
    
except Exception as e:
    print(f"[ERROR] Failed to write output: {e}", file=sys.stderr)
    sys.exit(1)
