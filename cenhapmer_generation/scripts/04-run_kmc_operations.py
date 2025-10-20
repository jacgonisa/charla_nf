#!/usr/bin/env python3

"""
Script: 04-run_kmc_operations.py
Description: Execute KMC complex operations to generate unique k-mer databases
             (cenhapmers) for each accession/chromosome/region

Usage: python3 04-run_kmc_operations.py [options]

Options:
    --input_dir DIR        Directory with .op files (default: 02-kmc_complex_op_files)
    --output_dir DIR       Output directory for cenhapmers (default: 03-cenhapmers)
    --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
    --parallel N           Number of parallel kmc_tools processes (default: 1)
    --help                 Show this help message

Requirements:
    - kmc_tools (from KMC package)
    - Operation files from step 03

Output:
    - ${output_dir}/k${K}/unique_${ACCESSION}_{CEN|ARMS}_${CHR}_k${K}.kmc_{pre,suf}

These are the final cenhapmer databases containing k-mers unique to each
accession/chromosome/region combination.
"""

import os
import argparse
import subprocess
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_kmc_tools(op_file, output_dir, verbose=True):
    """Run kmc_tools with an operation file and move results to output directory"""

    filename = Path(op_file).name
    # Remove .op extension and split by underscore
    parts = filename[:-3].split('_')

    # Extract info from filename: unique_ACCESSION_REGION_CHR_kSIZE.op
    if len(parts) < 5:
        print(f"⚠️  Warning: Unexpected filename format: {filename}")
        return False

    k_size = parts[-1]  # e.g., 'k21'
    accession = parts[1]
    region = parts[2]
    chrom = parts[3]

    if verbose:
        print(f"\n>>> Processing: {accession} {region} {chrom} ({k_size})")

    # Create output directory for this k-mer size
    kmc_dir = os.path.join(output_dir, k_size)
    os.makedirs(kmc_dir, exist_ok=True)

    # Base name for output database
    output_db_base = f"unique_{accession}_{region}_{chrom}_{k_size}"

    try:
        if verbose:
            print(f"    Running kmc_tools complex {op_file}")

        # Run kmc_tools
        result = subprocess.run(
            ["kmc_tools", "complex", op_file],
            capture_output=True,
            text=True,
            check=True
        )

        # Move output files to the appropriate directory
        kmc_pre_file = f"{output_db_base}.kmc_pre"
        kmc_suf_file = f"{output_db_base}.kmc_suf"

        if os.path.exists(kmc_pre_file) and os.path.exists(kmc_suf_file):
            dest_pre = os.path.join(kmc_dir, kmc_pre_file)
            dest_suf = os.path.join(kmc_dir, kmc_suf_file)

            shutil.move(kmc_pre_file, dest_pre)
            shutil.move(kmc_suf_file, dest_suf)

            if verbose:
                print(f"    ✅ Success! Output: {kmc_dir}/{output_db_base}.kmc_{{pre,suf}}")

            return True
        else:
            print(f"    ❌ Error: Output files not found for {output_db_base}")
            return False

    except subprocess.CalledProcessError as e:
        print(f"    ❌ Error running kmc_tools: {e}")
        if e.stderr:
            print(f"    Error output: {e.stderr}")
        return False
    except Exception as e:
        print(f"    ❌ Unexpected error: {e}")
        return False


def process_single_kmer_size(input_dir, output_dir, k_size, parallel, verbose=True):
    """Process all operation files for a single k-mer size"""

    k_dir = os.path.join(input_dir, f"k{k_size}")

    if not os.path.exists(k_dir):
        print(f"⚠️  Warning: Directory not found: {k_dir}")
        return 0, 0

    # Find all .op files
    op_files = []
    for root, dirs, files in os.walk(k_dir):
        for file in files:
            if file.endswith('.op'):
                op_files.append(os.path.join(root, file))

    if not op_files:
        print(f"⚠️  Warning: No .op files found in {k_dir}")
        return 0, 0

    print(f"\n{'='*60}")
    print(f"Processing k-mer size: {k_size}")
    print(f"Total operation files: {len(op_files)}")
    print(f"Parallel processes: {parallel}")
    print(f"{'='*60}")

    successful = 0
    failed = 0

    if parallel == 1:
        # Sequential processing
        for i, op_file in enumerate(op_files, 1):
            print(f"\n[{i}/{len(op_files)}] Processing {Path(op_file).name}")
            if run_kmc_tools(op_file, output_dir, verbose):
                successful += 1
            else:
                failed += 1
    else:
        # Parallel processing
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = {
                executor.submit(run_kmc_tools, op_file, output_dir, False): op_file
                for op_file in op_files
            }

            for i, future in enumerate(as_completed(futures), 1):
                op_file = futures[future]
                print(f"[{i}/{len(op_files)}] {Path(op_file).name}...", end=" ")
                try:
                    if future.result():
                        print("✅")
                        successful += 1
                    else:
                        print("❌")
                        failed += 1
                except Exception as e:
                    print(f"❌ Exception: {e}")
                    failed += 1

    return successful, failed


def main():
    parser = argparse.ArgumentParser(
        description="Run KMC complex operations to generate cenhapmers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--input_dir', default='02-kmc_complex_op_files',
                       help='Directory with .op files (default: 02-kmc_complex_op_files)')
    parser.add_argument('--output_dir', default='03-cenhapmers',
                       help='Output directory for cenhapmers (default: 03-cenhapmers)')
    parser.add_argument('--kmer_sizes', default='21,31,41',
                       help='Comma-separated k-mer sizes (default: 21,31,41)')
    parser.add_argument('--parallel', type=int, default=1,
                       help='Number of parallel kmc_tools processes (default: 1)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Parse k-mer sizes
    kmer_sizes = [int(k) for k in args.kmer_sizes.split(',')]

    print(f"\n{'='*60}")
    print(f"KMC Operations Execution")
    print(f"{'='*60}")
    print(f"Input directory:  {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"K-mer sizes:      {', '.join(map(str, kmer_sizes))}")
    print(f"Parallel jobs:    {args.parallel}")
    print(f"{'='*60}")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    total_successful = 0
    total_failed = 0

    # Process each k-mer size
    for k_size in kmer_sizes:
        successful, failed = process_single_kmer_size(
            args.input_dir,
            args.output_dir,
            k_size,
            args.parallel,
            args.verbose
        )
        total_successful += successful
        total_failed += failed

    # Final summary
    print(f"\n{'='*60}")
    print(f"FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Total successful: {total_successful}")
    print(f"Total failed:     {total_failed}")
    print(f"Success rate:     {total_successful/(total_successful+total_failed)*100:.1f}%")
    print(f"{'='*60}")
    print(f"\n✅ Cenhapmer generation complete!")
    print(f"Output directory: {args.output_dir}")

    if total_failed > 0:
        print(f"\n⚠️  Warning: {total_failed} operations failed")
        return 1
    return 0


if __name__ == "__main__":
    exit(main())
