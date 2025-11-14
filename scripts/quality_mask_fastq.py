#!/usr/bin/env python3
"""
Quality-based base masking for FASTQ files.
Replaces bases with quality scores below threshold with 'N'.
"""

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
import gzip

def replace_low_quality_bases(input_fastq, output_fastq, min_quality):
    """
    Reads a fastq file and replaces bases with a quality score lower than
    min_quality with 'N'.

    Args:
        input_fastq (str): The path to the input fastq file (plain or gzipped).
        output_fastq (str): The path to the output fastq file.
        min_quality (int): The minimum Phred quality score. Bases below this
                           score will be replaced by 'N'.
    """
    # Handle both gzipped and plain FASTQ files
    if input_fastq.endswith('.gz'):
        in_handle = gzip.open(input_fastq, "rt")
    else:
        in_handle = open(input_fastq, "r")

    with in_handle, open(output_fastq, "w") as out_handle:
        processed_reads = 0
        masked_bases = 0
        total_bases = 0
        
        for title, seq, qual in FastqGeneralIterator(in_handle):
            new_seq = ""
            new_qual = ""
            
            # Iterate through each base and its quality score
            for i in range(len(seq)):
                total_bases += 1
                # Convert ASCII quality score to integer Phred score
                quality_score = ord(qual[i]) - 33
                
                if quality_score < min_quality:
                    # Replace base with 'N' and keep the original quality score
                    new_seq += "N"
                    new_qual += qual[i]
                    masked_bases += 1
                else:
                    # Keep the original base and quality score
                    new_seq += seq[i]
                    new_qual += qual[i]
            
            # Write the modified record to the output file
            out_handle.write(f"@{title}\n{new_seq}\n+\n{new_qual}\n")
            processed_reads += 1
            
            if processed_reads % 10000 == 0:
                print(f"Processed {processed_reads} reads...", file=sys.stderr)
        
        # Calculate masking statistics
        masking_rate = (masked_bases / total_bases * 100) if total_bases > 0 else 0
        
        print(f"Quality masking complete:", file=sys.stderr)
        print(f"  Processed reads: {processed_reads:,}", file=sys.stderr)
        print(f"  Total bases: {total_bases:,}", file=sys.stderr)
        print(f"  Masked bases (Q<{min_quality}): {masked_bases:,}", file=sys.stderr)
        print(f"  Masking rate: {masking_rate:.2f}%", file=sys.stderr)

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 quality_mask_fastq.py <input.fastq> <output.fastq> <min_quality>")
        print("Example: python3 quality_mask_fastq.py input.fq output_masked.fq 20")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        min_qual = int(sys.argv[3])
    except ValueError:
        print("Error: Minimum quality must be an integer.")
        sys.exit(1)
    
    if min_qual < 0 or min_qual > 40:
        print(f"Warning: Unusual quality threshold ({min_qual}). Typical values are 10-30.")
    
    print(f"Replacing bases in {input_file} with Phred score < {min_qual} with 'N'...", file=sys.stderr)
    
    try:
        replace_low_quality_bases(input_file, output_file, min_qual)
        print(f"Quality-masked FASTQ file saved as: {output_file}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
