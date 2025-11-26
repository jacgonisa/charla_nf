import os
import sys
import subprocess

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <cpus>")
        print("Example: python script.py sample_nonhybrid 8")
        sys.exit(1)
    
    input_directory = sys.argv[1]
    cpus = sys.argv[2]
    
    # Extract sample_id from input directory (remove '_nonhybrid' suffix)
    sample_id = input_directory.replace('_nonhybrid', '')
    
    OUTPUT_ALIGNMENTS_DIR = os.path.join(input_directory, "mappings")

    REFERENCE_GENOMES = {
        "Col-0": "index/Col-0_renamed.fa",
        "Ler-0": "index/Ler-0_renamed.fa"
    }

    # Simplified to only 2 mapping modes for faster processing
    MODE_CMD = {
        "mm_ont": f"minimap2 --secondary=no -ax map-ont -t {cpus} -p 1.0 -N 5 ${{genome_ref}} ${{READS}} | samtools view -b -F 4 -q 10 - | samtools sort -o ${{OUTPUT_BAM}}",
        "wm_mapont": f"winnowmap -a --secondary=no -W ${{REPETITIVE_KMER_FILE}} -x map-ont -t {cpus} -p 1.0 -N 5 ${{genome_ref}} ${{READS}} | samtools view -b -F 4 -q 10 - | samtools sort -o ${{OUTPUT_BAM}}"
    }

    # Create mappings directory if it doesn't exist
    os.makedirs(OUTPUT_ALIGNMENTS_DIR, exist_ok=True)

    print(f"--- Starting Read Alignment ---")
    print(f"Input directory: {input_directory}/")
    print(f"Output alignments will be in: {OUTPUT_ALIGNMENTS_DIR}/")
    print(f"Using {cpus} CPU cores")

    # Check if input directory exists
    if not os.path.exists(input_directory):
        print(f"Error: Input directory '{input_directory}' not found!")
        sys.exit(1)

    # Check if reference genomes exist
    for genome_name, genome_path in REFERENCE_GENOMES.items():
        if not os.path.exists(genome_path):
            print(f"Warning: Reference genome '{genome_name}' not found at '{genome_path}'")

    fasta_files = [f for f in os.listdir(input_directory) if f.endswith("_sequences.fa")]
    
    if not fasta_files:
        print(f"No FASTA files found in {input_directory} ending with '_sequences.fa'")
        # Create empty BAM files to satisfy Nextflow output requirements
        for mode_name in MODE_CMD.keys():
            for haplotype in ['C', 'L']:
                empty_bam = os.path.join(OUTPUT_ALIGNMENTS_DIR, f"{haplotype}_{mode_name}.bam")
                with open(empty_bam, 'w') as f:
                    f.write("")
        print("Created empty BAM files for Nextflow output requirements")
        sys.exit(0)

    for fasta_file in fasta_files:
        haplotype_abbr = fasta_file.split('_')[0]
        full_fasta_path = os.path.join(input_directory, fasta_file)

        reference_genome_path = None
        genome_name = None
        if haplotype_abbr.startswith('C'):
            reference_genome_path = REFERENCE_GENOMES.get("Col-0")
            genome_name = "Col-0"
        elif haplotype_abbr.startswith('L'):
            reference_genome_path = REFERENCE_GENOMES.get("Ler-0")
            genome_name = "Ler-0"
        else:
            print(f"Warning: Unknown haplotype abbreviation '{haplotype_abbr}'. Skipping {fasta_file}.")
            continue

        if not reference_genome_path or not os.path.exists(reference_genome_path):
            print(f"Error: Reference genome not found for '{genome_name}'. Skipping {fasta_file}.")
            continue

        if not os.path.exists(full_fasta_path) or os.path.getsize(full_fasta_path) == 0:
            print(f"Skipping {fasta_file} as it's empty or doesn't exist.")
            continue

        print(f"\nProcessing {fasta_file} against {genome_name}...")

        for mode_name, cmd_template in MODE_CMD.items():
            output_bam_name = f"{haplotype_abbr}_{mode_name}.bam"
            output_bam_path = os.path.join(OUTPUT_ALIGNMENTS_DIR, output_bam_name)

            current_cmd_template = cmd_template

            if "winnowmap" in current_cmd_template:
                repetitive_kmer_file = f"index/{genome_name}_repetitive_k15.txt"
                current_cmd_template = current_cmd_template.replace("${REPETITIVE_KMER_FILE}", repetitive_kmer_file)

                if not os.path.exists(repetitive_kmer_file):
                    print(f"  Warning: Repetitive kmer file for '{genome_name}' not found at '{repetitive_kmer_file}'.")

            command = (current_cmd_template.replace("${genome_ref}", reference_genome_path)
                                          .replace("${READS}", full_fasta_path)
                                          .replace("${OUTPUT_BAM}", output_bam_path))

            print(f"  Running mode: {mode_name}")

            try:
                process = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
                print(f"  Successfully aligned using {mode_name}. Output: {output_bam_path}")
                
                # Create BAM index
                index_cmd = f"samtools index {output_bam_path}"
                subprocess.run(index_cmd, shell=True, check=True)
                print(f"  Created index: {output_bam_path}.bai")
                
            except subprocess.CalledProcessError as e:
                print(f"  Error running {mode_name} for {fasta_file}: {e}")
                if e.stderr:
                    print(f"  Stderr: {e.stderr[:200]}...")
            except Exception as e:
                print(f"  Unexpected error for {mode_name}: {e}")

    print("\n--- All alignments complete. ---")

if __name__ == "__main__":
    main()
