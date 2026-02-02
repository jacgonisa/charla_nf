#!/usr/bin/env python3

"""
Script: 03-generate_kmc_operations.py
Description: Generate KMC operation (.op) files to find unique k-mers for each
             accession/chromosome/region combination

Usage: python3 03-generate_kmc_operations.py [options]

Options:
    --input_dir DIR        Directory with k-mer databases (default: 01-all_kmers_per_Chr_CEN_and_accesion)
    --output_dir DIR       Output directory for .op files (default: 02-kmc_complex_op_files)
    --parent1 NAME         Name of parent 1 accession (default: Col-0)
    --parent2 NAME         Name of parent 2 accession (default: Ler-0)
    --num_chr N            Number of chromosomes (default: 5)
    --chr_prefix PREFIX    Chromosome prefix (default: Chr)
    --kmer_sizes SIZES     Comma-separated k-mer sizes (default: 21,31,41)
    --help                 Show this help message

Requirements:
    - Python 3.6+
    - K-mer databases from step 02

Output:
    - ${output_dir}/k${K}/unique_${ACCESSION}_{CEN|ARMS}_${CHR}_k${K}.op

Operation files use KMC set operations to find unique k-mers:
    unique_X = set_X - (set_1 + set_2 + ... + set_N)
where set_X is the target and all other sets are subtracted
"""

import os
import argparse
from pathlib import Path


def generate_operation_files(input_dir, output_dir, parent1, parent2, num_chr,
                            chr_prefix, kmer_sizes, skip_cen=False):
    """Generate KMC operation files for finding unique k-mers"""

    # Generate accessions and chromosomes lists
    accessions = [parent1, parent2]
    chromosomes = [f"{chr_prefix}{i}" for i in range(1, num_chr + 1)]
    regions = ["CHR"] if skip_cen else ["CEN", "ARMS"]

    print(f"\n=== KMC Operation File Generation ===")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Accessions:       {', '.join(accessions)}")
    print(f"Chromosomes:      {', '.join(chromosomes)}")
    print(f"Regions:          {', '.join(regions)}")
    print(f"K-mer sizes:      {', '.join(map(str, kmer_sizes))}")
    print(f"Skip centromeres: {skip_cen}")
    print(f"=====================================\n")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    total_files = 0

    for k in kmer_sizes:
        print(f"\n📊 Processing k-mer size: {k}")

        # Create dictionary to store set numbers and file paths
        sets = {}
        counter = 1

        # Build the sets dictionary
        for acc in accessions:
            for chrom in chromosomes:
                for reg in regions:
                    # Path to KMC database (without .kmc_pre/.kmc_suf extension)
                    file_path = f"{input_dir}/{acc}/{chrom}/{reg}/{acc}_{reg}_{chrom}_k{k}"
                    sets[counter] = file_path
                    counter += 1

        print(f"   Total sets: {len(sets)}")

        # Create output directory for this k-mer size
        k_dir = os.path.join(output_dir, f"k{k}")
        os.makedirs(k_dir, exist_ok=True)

        # Generate operation files for each combination
        for acc in accessions:
            for reg in regions:
                for chrom in chromosomes:
                    # Find the target set number
                    target_set = None
                    for set_num, path in sets.items():
                        if acc in path and reg in path and chrom in path:
                            target_set = set_num
                            break

                    if target_set is None:
                        print(f"   ⚠️  Warning: Could not find set for {acc} {reg} {chrom}")
                        continue

                    # Generate the KMC operation: unique = target - (all others)
                    other_sets = [num for num in range(1, len(sets) + 1) if num != target_set]
                    set_diff = f"set{target_set} - (" + " + ".join(f"set{i}" for i in other_sets) + ")"

                    # Output .op filename
                    filename = os.path.join(k_dir, f"unique_{acc}_{reg}_{chrom}_k{k}.op")

                    # Write operation file
                    with open(filename, 'w') as f:
                        # Write INPUT section
                        f.write("INPUT:\n")
                        for set_num, path in sets.items():
                            f.write(f"set{set_num} = {path}\n")

                        # Write OUTPUT section
                        f.write("\nOUTPUT:\n")
                        f.write(f"unique_{acc}_{reg}_{chrom}_k{k} = {set_diff}\n")

                    total_files += 1

        print(f"   ✅ Generated {len(accessions) * len(regions) * len(chromosomes)} operation files for k={k}")

    print(f"\n✅ Complete! Generated {total_files} operation files total")
    print(f"Output directory: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate KMC operation files for finding unique k-mers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--input_dir', default='01-all_kmers_per_Chr_CEN_and_accesion',
                       help='Directory with k-mer databases (default: 01-all_kmers_per_Chr_CEN_and_accesion)')
    parser.add_argument('--output_dir', default='02-kmc_complex_op_files',
                       help='Output directory for .op files (default: 02-kmc_complex_op_files)')
    parser.add_argument('--parent1', default='Col-0',
                       help='Name of parent 1 accession (default: Col-0)')
    parser.add_argument('--parent2', default='Ler-0',
                       help='Name of parent 2 accession (default: Ler-0)')
    parser.add_argument('--num_chr', type=int, default=5,
                       help='Number of chromosomes (default: 5)')
    parser.add_argument('--chr_prefix', default='Chr',
                       help='Chromosome prefix (default: Chr)')
    parser.add_argument('--kmer_sizes', default='21,31,41',
                       help='Comma-separated k-mer sizes (default: 21,31,41)')
    parser.add_argument('--skip-cen', action='store_true',
                       help='Skip centromere annotation (generate chromosome-level hapmers only)')

    args = parser.parse_args()

    # Parse k-mer sizes
    kmer_sizes = [int(k) for k in args.kmer_sizes.split(',')]

    generate_operation_files(
        args.input_dir,
        args.output_dir,
        args.parent1,
        args.parent2,
        args.num_chr,
        args.chr_prefix,
        kmer_sizes,
        skip_cen=args.skip_cen
    )


if __name__ == "__main__":
    main()
