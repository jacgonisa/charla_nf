import csv
import re
import sys

def create_curated_profile(raw, readmer_cutoff):
    """
    Processes a raw profile string and returns a curated profile string.
    Entries with quality score lower than 11 or with abbreviation '0' are removed.
    """
    entries = raw.split(',')
    curated = []
    for entry in entries:
        try:
            code, quality = entry.split(':')
            quality = int(quality)
        except ValueError:
            continue
        if quality < readmer_cutoff:
            continue
        if code == "0":
            continue
        curated.append(f"{code}:{quality}")
    return ",".join(curated)

def create_supercurated_profile(curated_str, threshold):
    """
    Merges consecutive entries with the same code that meet the threshold.
    """
    if not curated_str:
        return ""
    
    curated = [tuple(x.split(':')) for x in curated_str.split(',')]
    supercurated_groups = []
    current_code = curated[0][0]
    count = 1

    for entry in curated[1:]:
        code, _ = entry
        if code == current_code:
            count += 1
        else:
            if count >= threshold:
                supercurated_groups.append(f"{count}{current_code}")
            current_code = code
            count = 1

    if count >= threshold:
        supercurated_groups.append(f"{count}{current_code}")

    return ",".join(supercurated_groups)

def hybrid_label(supercurated_profile):
    """
    Creates the hybrid label based on the number of distinct codes.
    If 2 or more distinct codes exist, the label is "yes|<number>".
    Otherwise, it is "no|<number>".
    """
    if not supercurated_profile:
        return "no|0"
    
    groups = supercurated_profile.split(',')
    distinct_codes = set()
    
    for group in groups:
        m = re.match(r'(\d+)([A-Z0-9]+)', group)
        if m:
            _, code = m.groups()
            distinct_codes.add(code)
    
    count_codes = len(distinct_codes)
    if count_codes >= 2:
        return f"yes|{count_codes}"
    else:
        return f"no|{count_codes}"


def create_ultracurated_profile(supercurated_str):
    """
    Processes the supercurated profile into an ultra-curated profile by combining consecutive entries 
    with the same code into one, summing their counts.
    """
    if not supercurated_str:
        print("Supercurated string is empty.")
        return ""
    
    curated = []
    for entry in supercurated_str.split(','):
        # Match the count (numeric part) followed by the code (alphanumeric part)
        m = re.match(r'(\d+)([A-Z0-9]+)', entry)
        if m:
            count, code = m.groups()
            curated.append((int(count), code))
        else:
            # If the format doesn't match, skip this entry
            print(f"Skipping invalid entry: {entry}")
            continue
    
    if not curated:
        print("Curated list is empty.")
        return ""
    
    # Now process the curated list to create the ultra-curated profile
    ultra_curated_groups = []
    current_code = curated[0][1]
    total_count = curated[0][0]

    for entry in curated[1:]:
        count, code = entry
        if code == current_code:
            total_count += count
        else:
            ultra_curated_groups.append(f"{total_count}{current_code}")
            current_code = code
            total_count = count

    # Add the last group
    ultra_curated_groups.append(f"{total_count}{current_code}")
    
    return ",".join(ultra_curated_groups)


#import re


import re

def classify_hybrid(ultracurated_profile):
    """
    Further classifies hybrid reads (based on the ultra-curated profile, which is assumed
    to contain only entries with a count followed immediately by a code, separated by commas)
    into:
      - SCO (Simple CrossOver): exactly 2 groups (e.g., "35LC1,15CC1")
      - NCO (Non CrossOver): exactly 3 groups in an ABA pattern (e.g., "35LC1,15CC1,40LC1")
      - GC-CO (Gene Conversion-associated CO): more than 3 groups that strictly alternate 
        between the two codes (e.g., "35LC1,5CC1,6LC1,27CC1")
    Any sequence that does not meet these criteria (or does not have exactly 2 distinct codes) 
    is returned as "Unknown".
    """
    if not ultracurated_profile:
        return "None"

    # Split into groups and parse each group into (count, code)
    groups = ultracurated_profile.split(',')
    parsed = []
    for g in groups:
        m = re.match(r'(\d+)([A-Z0-9]+)', g)
        if m:
            count, code = m.groups()
            parsed.append((int(count), code))
    
    if not parsed:
        return "None"
    
    # Extract the sequence of codes
    codes = [code for _, code in parsed]
    
    # Must have exactly 2 distinct codes to be considered
    if len(set(codes)) != 2:
        return "Unknown"
    
    n = len(parsed)
    # If exactly 2 groups, classify as SCO.
    if n == 2:
        return "SCO"
    
    # Check that the sequence is alternating: each adjacent group must have a different code.
    alternating = all(codes[i] != codes[i+1] for i in range(n-1))
    if not alternating:
        return "Unknown"
    
    # If there are exactly 3 groups and the pattern is ABA, classify as NCO.
    if n == 3:
        if codes[0] == codes[2]:
            return "NCO"
        else:
            return "Unknown"
    
    # If there are 4 or more groups and they strictly alternate, classify as GC-CO.
    if n == 4:
        return "GC-CO"
    
    return "Unknown"



def ectopic_classification(code1, code2):
    """
    Classifies ectopic recombination between two codes (for yes|2 hybrid reads).
    Codes are expected to be in the format two letters followed by digits (e.g., LC1, CA2).
    
    The logic is as follows:
      1. If the substring after the first letter is identical (i.e. code1[1:] == code2[1:]),
         they are homologous → return "no".
      2. If they come from the same accession (same first letter) and map to the same chromosome 
         (i.e. numeric parts are equal) but the two-letter prefixes differ, return "no|edge".
      3. Otherwise, the recombination is ectopic, and we further classify as:
             yes|intraaccession|<chromosome_status>   or
             yes|interaccession|<chromosome_status>
         where <chromosome_status> is "intrachromosome" if the numeric parts are equal, and
         "interchromosome" if they differ.
    """
    m1 = re.match(r'([A-Z]{2})(\d+)', code1)
    m2 = re.match(r'([A-Z]{2})(\d+)', code2)
    if not m1 or not m2:
        return "None"
    
    prefix1, n1 = m1.groups()
    prefix2, n2 = m2.groups()
    num1, num2 = int(n1), int(n2)
    
    # Homologous if they differ only by the first letter.
    if code1[1:] == code2[1:]:
        return "no"
    
    # Edge case: same accession and same chromosome, but different two-letter prefixes.
    if code1[0] == code2[0] and num1 == num2 and prefix1 != prefix2:
        return "no|edge"
    
    # Ectopic:
    access = "intraaccession" if code1[0] == code2[0] else "interaccession"
    chrom = "intrachromosome" if num1 == num2 else "interchromosome"
    return f"yes|{access}|{chrom}"

def process_read(readname, raw_profile, thresholds, readmer_cutoff):
    """
    Given a read name and its raw profile, processes it across multiple thresholds.
    Returns a list of tuples with:
      (readname, threshold, curated_profile, super_curated_profile, ultra_curated_profile, hybrid, hybrid_type, ectopic)
    """
    curated = create_curated_profile(raw_profile, readmer_cutoff)
    results = []
    for threshold in thresholds:
        supercurated = create_supercurated_profile(curated, threshold)
        ultra_curated = create_ultracurated_profile(supercurated)  # Generate ultra curated profile
        hybrid = hybrid_label(ultra_curated)  # Use ultra curated for hybrid classification
        
        # Further classify hybrid types only if exactly 2 distinct codes in ultra curated profile
        if hybrid.startswith("yes|2"):
            hybrid_type = classify_hybrid(ultra_curated)  # Classify based on ultra curated
            # Extract distinct codes from the ultra curated profile:
            codes_set = set()
            for group in ultra_curated.split(','):
                m = re.match(r'\d+([A-Z0-9]+)', group)
                if m:
                    codes_set.add(m.group(1))
            if len(codes_set) == 2:
                code1, code2 = sorted(list(codes_set))
                ectopic = ectopic_classification(code1, code2)
            else:
                ectopic = "None"
        else:
            hybrid_type = "None"
            ectopic = "None"
        
        results.append((readname, threshold, curated, supercurated, ultra_curated, hybrid, hybrid_type, ectopic))
    return results

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    readmer_cutoff = int(sys.argv[3].strip())
    thresholds = [10,20,30,35,40,45,50,75,100]

    # Open the output file for writing before processing any data
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["readname", "threshold", "curated_profile", "super_curated_profile", "ultra_curated_profile", "hybrid", "hybrid_type", "ectopic"])
        
        with open(input_file, "r") as infile:
            lines = infile.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith(">"):
                readname = line[1:].strip()
                if i + 1 < len(lines):
                    raw_profile = lines[i + 1].strip()
                    # Process the read and write the results progressively
                    for row in process_read(readname, raw_profile, thresholds, readmer_cutoff):
                        writer.writerow(row)
                i += 2
            else:
                i += 1

    print(f"TSV file '{output_file}' created.")
if __name__ == "__main__":
    main()
