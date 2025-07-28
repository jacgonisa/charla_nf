#!/usr/bin/env python3
import sys

def process_kmer_profile(read_name, profile, k=81):
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

def merge_consecutive_by_code(segments):
    """
    Given a list of segments (each a tuple:
      (read, start, end, region, token_count)) in the original token order,
    first filter out any segment that is "isolated and small" (i.e. length < 100 bp 
    and fewer than 5 consecutive tokens) and then merge all segments that are consecutive 
    (in the filtered order) and have the same region code.
    
    Merging simply takes the start of the first segment and the end of the last segment.
    """
    if not segments:
        return []
    
    # Filter out isolated small segments.
    filtered = [
        seg for seg in segments
        if not ((seg[2] - seg[1]) < 100 and seg[4] < 100)
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

def main(kmer_profile_file, curated_tsv_file, outfile, k=81):
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
                    merged = merge_consecutive_by_code(segments)
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
            merged = merge_consecutive_by_code(segments)
            for rec in merged:
                r, start, end, region, _ = rec
                bed_lines.append(f"{r}\t{start}\t{end}\t{region}")
    
    with open(outfile, "w") as out:
        for line in bed_lines:
            out.write(line + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 script.py <input_kmer_profile.fa> <curated_tsv_file> <output.bed>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])


