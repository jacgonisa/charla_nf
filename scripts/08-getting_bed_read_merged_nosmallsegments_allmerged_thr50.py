#!/usr/bin/env python3
"""
Generate BED file from k-mer profiles with configurable filtering thresholds.

This script processes k-mer profiles to create BED coordinates for read segments,
filtering out small segments based on configurable length and token count thresholds.
"""
import sys
import argparse

def process_kmer_profile(read_name, profile, k):
    """
    Given a read name and its kmer profile string, group contiguous tokens 
    (ignoring the numeric value after the colon) and return a list of tuples:
      (read_name, block_start, block_end, region, token_count)
    Each token covers positions [i, i+k) and a block spanning tokens i to j 
    is represented as [i, j+k). The token_count is the number of consecutive tokens.
    """
    tokens = [token.strip() for token in profile.split(',') if token.strip()]
    if not tokens:
        return []
    
    # Extract only the region label from each token (ignoring any number after the colon)
    labels = [token.split(':')[0] for token in tokens]
    
    segments = []
    current_label = labels[0]
    start_index = 0

    for i, label in enumerate(labels):
        if label != current_label:
            block_start = start_index
            block_end = (i - 1) + k  # token (i-1) covers bases [i-1, i-1+k)
            token_count = i - start_index
            segments.append((read_name, block_start, block_end, current_label, token_count))
            current_label = label
            start_index = i
    # Append the final block.
    block_start = start_index
    block_end = (len(labels) - 1) + k
    token_count = len(labels) - start_index
    segments.append((read_name, block_start, block_end, current_label, token_count))
    
    return segments

def merge_consecutive_by_code(segments, min_segment_length=50, min_token_count=50):
    """
    Given a list of segments (each a tuple:
      (read, start, end, region, token_count)) in the original token order,
    first filter out any segment that is "isolated and small" and then merge all
    segments that are consecutive (in the filtered order) and have the same region code.

    A segment is filtered if BOTH conditions are true:
      - Length < min_segment_length bp AND
      - Token count < min_token_count

    Merging simply takes the start of the first segment and the end of the last segment.

    Parameters:
    -----------
    segments : list
        List of segment tuples (read, start, end, region, token_count)
    min_segment_length : int
        Minimum segment length in bp (default: 50)
    min_token_count : int
        Minimum number of consecutive k-mer tokens (default: 50)

    Returns:
    --------
    list : Filtered and merged segments
    """
    if not segments:
        return []

    # Filter out isolated small segments.
    # A segment is removed only if BOTH length AND token count are below thresholds
    filtered = [
        seg for seg in segments
        if not ((seg[2] - seg[1]) < min_segment_length and seg[4] < min_token_count)
    ]
    if not filtered:
        return []
    
    merged = []
    # Start with the first segment in filtered order.
    current = filtered[0]
    
    for seg in filtered[1:]:
        # If the new segment has the same region code, merge it with the current chain.
        if seg[3] == current[3]:
            # Merge: keep the original start and update the end to the new segment's end.
            # (Token counts are summed, though they are not used in the final output.)
            current = (current[0], current[1], seg[2], current[3], current[4] + seg[4])
        else:
            merged.append(current)
            current = seg
    merged.append(current)
    return merged

def read_curated_tsv(tsv_file):
    """
    Reads a TSV file where:
      - The first column is the read name.
      - The last column is a dash-separated list of allowed kmer codes (e.g. "CC1-LC1").
    Returns a dictionary mapping read_name -> set(allowed_codes).
    """
    curated = {}
    with open(tsv_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            read_name = parts[0]
            allowed_codes = parts[-1].split('-')
            curated[read_name] = set(allowed_codes)
    return curated

def main(kmer_profile_file, curated_tsv_file, outfile, k, min_segment_length=50, min_token_count=50):
    """
    Main processing function.

    Parameters:
    -----------
    kmer_profile_file : str
        Path to k-mer profile FASTA file
    curated_tsv_file : str
        Path to curated TSV file with allowed codes per read
    outfile : str
        Path to output BED file
    k : int
        K-mer size
    min_segment_length : int
        Minimum segment length in bp (default: 50)
    min_token_count : int
        Minimum number of consecutive k-mer tokens (default: 50)
    """
    curated = read_curated_tsv(curated_tsv_file)

    with open(kmer_profile_file) as f:
        lines = f.read().splitlines()

    bed_lines = []
    current_read = None
    profile_lines = []

    for line in lines:
        if line.startswith(">"):
            # Process the previous read, if any.
            if current_read is not None and profile_lines:
                profile = "".join(profile_lines)
                if current_read in curated:
                    allowed = curated[current_read]
                    segments = process_kmer_profile(current_read, profile, k)
                    # Keep only segments with allowed codes.
                    segments = [seg for seg in segments if seg[3] in allowed]
                    merged = merge_consecutive_by_code(segments, min_segment_length, min_token_count)
                    for rec in merged:
                        r, start, end, region, _ = rec
                        bed_lines.append(f"{r}\t{start}\t{end}\t{region}")
            current_read = line[1:].strip()
            profile_lines = []
        else:
            profile_lines.append(line.strip())

    # Process the final read.
    if current_read is not None and profile_lines:
        profile = "".join(profile_lines)
        if current_read in curated:
            allowed = curated[current_read]
            segments = process_kmer_profile(current_read, profile, k)
            segments = [seg for seg in segments if seg[3] in allowed]
            merged = merge_consecutive_by_code(segments, min_segment_length, min_token_count)
            for rec in merged:
                r, start, end, region, _ = rec
                bed_lines.append(f"{r}\t{start}\t{end}\t{region}")

    with open(outfile, "w") as out:
        for line in bed_lines:
            out.write(line + "\n")

    print(f"Processed {len(bed_lines)} segments with thresholds: "
          f"min_length={min_segment_length}bp, min_tokens={min_token_count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate BED file from k-mer profiles with configurable filtering thresholds',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Use default thresholds (50bp, 50 tokens)
    python3 script.py profile.fa curated.tsv output.bed 21

    # Use custom thresholds (25bp, 25 tokens for more sensitive detection)
    python3 script.py profile.fa curated.tsv output.bed 21 --min-length 25 --min-tokens 25

    # Use very strict filtering (100bp, 100 tokens)
    python3 script.py profile.fa curated.tsv output.bed 21 --min-length 100 --min-tokens 100
        """
    )

    parser.add_argument('kmer_profile', help='Input k-mer profile FASTA file')
    parser.add_argument('curated_tsv', help='Curated TSV file with allowed codes per read')
    parser.add_argument('output_bed', help='Output BED file')
    parser.add_argument('k_size', type=int, help='K-mer size')
    parser.add_argument('--min-length', type=int, default=50,
                       help='Minimum segment length in bp (default: 50)')
    parser.add_argument('--min-tokens', type=int, default=50,
                       help='Minimum number of consecutive k-mer tokens (default: 50)')

    args = parser.parse_args()

    main(args.kmer_profile, args.curated_tsv, args.output_bed, args.k_size,
         args.min_length, args.min_tokens)

