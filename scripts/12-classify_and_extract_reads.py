import os
import sys
import subprocess
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Classify and extract non-hybrid reads')
    parser.add_argument('kmer_matrix_file', help='TSV file with kmer counts per read')
    parser.add_argument('main_fasta_file', help='Main FASTA file with all reads')
    parser.add_argument('output_base_dir', help='Output directory')
    parser.add_argument('kmer_size', type=int, help='Kmer size value')
    
    args = parser.parse_args()
    
    # Use the parsed arguments
    input_kmer_filepath = args.kmer_matrix_file
    main_fasta_filepath = args.main_fasta_file
    OUTPUT_BASE_DIR = args.output_base_dir
    kmer_size = args.kmer_size

    # Create output directory
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

    # Mapping from full column names to abbreviated codes - WITH .txt extension
    col_to_abbr = {
        f"Col_ARMS_Chr1_k{kmer_size}.txt": 'CA1', f'Col_ARMS_Chr2_k{kmer_size}.txt': 'CA2', f'Col_ARMS_Chr3_k{kmer_size}.txt': 'CA3',
        f'Col_ARMS_Chr4_k{kmer_size}.txt': 'CA4', f'Col_ARMS_Chr5_k{kmer_size}.txt': 'CA5',
        f'Col_CEN_Chr1_k{kmer_size}.txt': 'CC1', f'Col_CEN_Chr2_k{kmer_size}.txt': 'CC2', f'Col_CEN_Chr3_k{kmer_size}.txt': 'CC3',
        f'Col_CEN_Chr4_k{kmer_size}.txt': 'CC4', f'Col_CEN_Chr5_k{kmer_size}.txt': 'CC5',
        f'Ler_ARMS_Chr1_k{kmer_size}.txt': 'LA1', f'Ler_ARMS_Chr2_k{kmer_size}.txt': 'LA2', f'Ler_ARMS_Chr3_k{kmer_size}.txt': 'LA3',
        f'Ler_ARMS_Chr4_k{kmer_size}.txt': 'LA4', f'Ler_ARMS_Chr5_k{kmer_size}.txt': 'LA5',
        f'Ler_CEN_Chr1_k{kmer_size}.txt': 'LC1', f'Ler_CEN_Chr2_k{kmer_size}.txt': 'LC2', f'Ler_CEN_Chr3_k{kmer_size}.txt': 'LC3',
        f'Ler_CEN_Chr4_k{kmer_size}.txt': 'LC4', f'Ler_CEN_Chr5_k{kmer_size}.txt': 'LC5'
    }

    # Rest of the classification logic stays the same
    haplotype_reads = {abbr: [] for abbr in col_to_abbr.values()}
    print(f"--- Starting Read Classification from: {input_kmer_filepath} ---")

    try:
        with open(input_kmer_filepath, 'r') as f:
            lines = f.readlines()

        header = lines[0].split()
        data_rows = lines[1:]
        kmer_col_indices = {header[i]: i for i in range(1, len(header))}

        for row in data_rows:
            parts = row.split()
            read_name = parts[0]
            current_read_kmer_counts = {}

            for col_name, abbr in col_to_abbr.items():
                if col_name in kmer_col_indices:
                    current_read_kmer_counts[abbr] = int(parts[kmer_col_indices[col_name]])

            dominant_haplotype_abbr = None
            haplotypes_over_50 = []

            for abbr, count in current_read_kmer_counts.items():
                if count > 50:
                    haplotypes_over_50.append(abbr)

            if len(haplotypes_over_50) == 1:
                potential_dominant = haplotypes_over_50[0]
                other_haplotype_kmer_sum = 0

                for abbr, count in current_read_kmer_counts.items():
                    if abbr != potential_dominant:
                        other_haplotype_kmer_sum += count

                if other_haplotype_kmer_sum <= 10:
                    dominant_haplotype_abbr = potential_dominant

            if dominant_haplotype_abbr:
                haplotype_reads[dominant_haplotype_abbr].append(read_name)

        print(f"Read classification complete. Writing read name lists to '{OUTPUT_BASE_DIR}/'.")

        for abbr, reads in haplotype_reads.items():
            read_list_filepath = os.path.join(OUTPUT_BASE_DIR, f"{abbr}_reads.txt")
            with open(read_list_filepath, 'w') as f:
                for read_name in reads:
                    f.write(f"{read_name}\n")

        print(f"\n--- Starting Sequence Extraction from: {main_fasta_filepath} ---")
        print(f"Output sequences will be saved to: {OUTPUT_BASE_DIR}/\n")

        for abbr in col_to_abbr.values():
            read_list_file = os.path.join(OUTPUT_BASE_DIR, f"{abbr}_reads.txt")
            output_fasta_file = os.path.join(OUTPUT_BASE_DIR, f"{abbr}_sequences.fa")

            if os.path.exists(read_list_file) and os.path.getsize(read_list_file) > 0:
                print(f"Extracting sequences for {abbr}...")
                result = subprocess.run(
                    ['seqkit', 'grep', '-f', read_list_file, main_fasta_filepath, '-o', output_fasta_file],
                    capture_output=True, text=True, check=False
                )
                if result.returncode == 0:
                    print(f"  Successfully created {output_fasta_file}")
                else:
                    print(f"  Error extracting sequences for {abbr}: {result.stderr}")
            else:
                print(f"  No reads found for {abbr}. Skipping sequence extraction.")

        print("\n--- Process Complete ---")

    except FileNotFoundError:
        print(f"Error: One of the input files not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
