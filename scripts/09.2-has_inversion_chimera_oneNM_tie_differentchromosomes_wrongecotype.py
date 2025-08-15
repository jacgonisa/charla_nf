#!/usr/bin/env python3
import sys
import os
import argparse

def parse_read_id_structure(read_id):
    """
    Parse a read ID to extract base_id and segment info.
    Expected format: 'base_id_segmentIndex_segmentLabel'
    
    The challenge is that base_id itself may contain underscores.
    We assume the last two underscore-separated parts are segmentIndex and segmentLabel,
    and everything before that is the base_id.
    
    Examples:
    - 'read_981729_len6991_err1.00pc_orig_start115774872_2_LA5' 
      -> base_id: 'read_981729_len6991_err1.00pc_orig_start115774872', segment: '2_LA5'
    - 'f389b71d-80f2-4cbd-98e6-2fd3a413ce69_3_CA2'
      -> base_id: 'f389b71d-80f2-4cbd-98e6-2fd3a413ce69', segment: '3_CA2'
    - 'simple_read_full' (no segment info)
      -> base_id: 'simple_read_full', segment: None
    """
    parts = read_id.split("_")
    
    if len(parts) < 3:
        # No segment info, entire read_id is base_id
        return read_id, None
    
    # Check if the last two parts look like segment info (index + label)
    # The segment index should be numeric, and label should start with L or C
    potential_index = parts[-2]
    potential_label = parts[-1]
    
    # Simple heuristic: if second-to-last part is numeric and last part starts with L or C
    try:
        int(potential_index)
        if potential_label.startswith(('L', 'C')):
            # This looks like segment info
            base_id = "_".join(parts[:-2])
            segment = f"{potential_index}_{potential_label}"
            return base_id, segment
    except ValueError:
        pass
    
    # If we can't identify segment info, treat entire read_id as base_id
    return read_id, None

def parse_paf_line(line):
    """
    Parse a PAF line and extract key fields.
    Also compute a 'base_id' using the improved parsing.
    """
    fields = line.strip().split("\t")
    if len(fields) < 12:
        return None
    rec = {}
    rec['read_id'] = fields[0]
    rec['base_id'], segment_info = parse_read_id_structure(rec['read_id'])
    
    try:
        rec['qlen'] = int(fields[1])
        rec['qstart'] = int(fields[2])
        rec['qend'] = int(fields[3])
    except ValueError:
        rec['qlen'] = rec['qstart'] = rec['qend'] = 0
    rec['strand'] = fields[4]
    rec['target'] = fields[5]
    try:
        rec['tlen'] = int(fields[6])
        rec['tstart'] = int(fields[7])
        rec['tend'] = int(fields[8])
    except ValueError:
        rec['tlen'] = rec['tstart'] = rec['tend'] = 0
    try:
        rec['match'] = int(fields[9])
        rec['aln_block'] = int(fields[10])
        rec['mapq'] = int(fields[11])
    except ValueError:
        rec['match'] = rec['aln_block'] = rec['mapq'] = 0
    for tag in fields[12:]:
        if tag.startswith("NM:i:"):
            try:
                rec['NM'] = int(tag[5:])
            except ValueError:
                rec['NM'] = None
        elif tag.startswith("AS:i:"):
            try:
                rec['AS'] = int(tag[5:])
            except ValueError:
                rec['AS'] = None
        elif tag.startswith("rl:i:"):
            try:
                rec['rl'] = int(tag[5:])
            except ValueError:
                rec['rl'] = None
    if 'AS' not in rec or 'NM' not in rec:
        return None
    rec['aln_qlen'] = rec['qend'] - rec['qstart']
    rec['line'] = line
    return rec

def get_segment_label(read_id):
    """
    Extract the segment label from a read_id using the improved parsing.
    Returns the segment part (e.g., '2_LA5') or None if no segment info.
    """
    base_id, segment_info = parse_read_id_structure(read_id)
    return segment_info

def get_read_ecotype(read_id):
    """
    Determines the ecotype (Col or Ler) from the segment label in read_id.
    Uses the improved parsing to extract segment info.
    Returns 'Col', 'Ler', or None if not clearly identifiable.
    """
    base_id, segment_info = parse_read_id_structure(read_id)
    
    if segment_info is None:
        return None
        
    # Extract the label part (after the last underscore in segment_info)
    segment_parts = segment_info.split("_")
    if len(segment_parts) >= 2:
        segment_label_suffix = segment_parts[1]  # e.g., 'LA5', 'CA2'
        if segment_label_suffix.startswith('C'):  # e.g., CA1, CA2
            return "Col"
        elif segment_label_suffix.startswith('L'):  # e.g., LA1, LA2
            return "Ler"
        elif "Col" in segment_label_suffix:  # For cases where the label directly says 'Col'
            return "Col"
        elif "Ler" in segment_label_suffix:  # For cases where the label directly says 'Ler'
            return "Ler"
    return None

# Test function to verify the parsing works correctly
def test_read_id_parsing():
    """Test the read ID parsing with various formats."""
    test_cases = [
        "read_981729_len6991_err1.00pc_orig_start115774872_2_LA5",
        "f389b71d-80f2-4cbd-98e6-2fd3a413ce69_3_CA2",
        "simple_read_1_Col",
        "no_segment_info",
        "another_complex_read_name_with_underscores_5_LB3"
    ]
    
    print("Testing read ID parsing:")
    for read_id in test_cases:
        base_id, segment = parse_read_id_structure(read_id)
        ecotype = get_read_ecotype(read_id)
        print(f"  {read_id}")
        print(f"    -> base_id: '{base_id}'")
        print(f"    -> segment: '{segment}'")
        print(f"    -> ecotype: '{ecotype}'")
        print()
   

def load_alignments(file_paths):
    """
    Read all alignments from the provided PAF files and group them by base read id.
    """
    alignments = {}
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            for line in f:
                rec = parse_paf_line(line)
                if rec is None:
                    continue
                base_id = rec['base_id']
                alignments.setdefault(base_id, []).append(rec)
    print(f"Loaded alignments for {len(alignments)} base read IDs.")
    return alignments

def get_inversion_details(aln_list):
    """
    Identify inversion evidence among segments.
    Returns a details string listing the inversion pair(s) found.
    """
    details = []
    n = len(aln_list)
    for i in range(n):
        for j in range(i+1, n):
            if aln_list[i]['target'] != aln_list[j]['target']:
                continue
            if aln_list[i]['strand'] == aln_list[j]['strand']:
                continue
            overlap_start = max(aln_list[i]['tstart'], aln_list[j]['tstart'])
            overlap_end = min(aln_list[i]['tend'], aln_list[j]['tend'])
            if overlap_start < overlap_end:
                details.append(f"(Segment {aln_list[i]['read_id']} vs {aln_list[j]['read_id']}: overlap {overlap_start}-{overlap_end})")
    return "; ".join(details)

def has_inversion(aln_list):
    """Return True if inversion evidence exists."""
    return bool(get_inversion_details(aln_list))

def similarity_score(rec):
    """
    Compute similarity based on the metric AS - NM.
    """
    if rec['AS'] is not None and rec['NM'] is not None:
        return rec['AS'] - rec['NM']
    return 0

def select_best_alignment(aln_list):
    """
    Select the best alignment from a list based on AS - NM.
    In case of a tie, choose the one with a higher fraction of the query aligned.
    This function does NOT consider expected ecotype for selection, only raw alignment quality.
    """
    best = None
    for rec in aln_list:
        score = similarity_score(rec)
        fraction = (rec['qend'] - rec['qstart']) / rec['qlen'] if rec['qlen'] > 0 else 0
        
        if best is None:
            best = rec
            best_score = score
            best_fraction = fraction
        else:
            if score > best_score or (score == best_score and fraction > best_fraction):
                best = rec
                best_score = score
                best_fraction = fraction
    return best

def group_by_target(aln_list):
    """
    Group alignments by target ecotype.
    If the target contains 'Col', assign group 'Col';
    if it contains 'Ler', assign group 'Ler'; otherwise 'Other'.
    """
    groups = {}
    for rec in aln_list:
        tgt = rec['target']
        if "Col" in tgt:
            grp = "Col"
        elif "Ler" in tgt:
            grp = "Ler"
        else:
            grp = "Other"
        groups.setdefault(grp, []).append(rec)
    return groups

def parse_target_arm(target):
    """
    Parse the target string to extract chromosome and arm.
    Expected format: 'ChrX_ARMS_..._Y' for arms.
    Returns (chromosome, arm) if possible; otherwise (target, None).
    """
    parts = target.split("_")
    if len(parts) >= 4 and "ARMS" in target:
        chrom = parts[0]
        arm = parts[-1]
        return chrom, arm
    else:
        chrom = parts[0]
        return chrom, None

def is_different_chromosomes_chimera(aln_list):
    """
    Check if a read's alignments map to different chromosomes.
    Returns (True, details) if more than one distinct chromosome is found.
    """
    chrom_set = set()
    details = []
    for rec in aln_list:
        chrom = rec['target'].split("_")[0]
        chrom_set.add(chrom)
        details.append(f"(chromosome: {chrom})")
    if len(chrom_set) > 1:
        return (True, "Different chromosomes: " + ", ".join(sorted(chrom_set)))
    return (False, "Same chromosome")

def is_different_arms_chimera(aln_list):
    """
    Check if a read's alignments map to different arms of the same chromosome.
    Returns (True, details) if, for any chromosome, alignments map to different arms.
    """
    chrom_arms = {}
    details = []
    for rec in aln_list:
        chrom, arm = parse_target_arm(rec['target'])
        if arm is not None:
            chrom_arms.setdefault(chrom, set()).add(arm)
    chim_flag = False
    for chrom, arms in chrom_arms.items():
        if len(arms) > 1:
            chim_flag = True
            details.append(f"(Chromosome {chrom}: arms {', '.join(sorted(arms))})")
    return (chim_flag, "; ".join(details))

def is_different_chromosomes(best_matches):
    """
    Check if the best alignments (one per segment group) map to different chromosomes.
    best_matches is a list of records.
    Returns (True, details) if more than one distinct chromosome is found.
    """
    chrom_set = set()
    details = []
    for rec in best_matches:
        chrom = rec['target'].split("_")[0]
        chrom_set.add(chrom)
        details.append(f"(chromosome: {chrom})")
    if len(chrom_set) > 1:
        return (True, "Different chromosomes: " + ", ".join(sorted(chrom_set)))
    return (False, "Same chromosome")

def process_read_segments(aln_list, coverage_tolerance=0.1):
    """
    For a given read, group alignments by segment label (using the segment index and label).
    For each segment group:
        - Identify the best alignment based on raw score (AS - NM).
        - If multiple alignments exist, check for one-difference or inconclusive ties.
    Returns:
        - best_segs: dictionary of best alignment(s) per segment group.
        - other_segs: dictionary of remaining alignments per segment group.
        - one_diff_segs: dictionary of segment groups with one-difference condition.
        - tie_segs: dictionary of segment groups with inconclusive ties.
    """
    seg_groups = {}
    for rec in aln_list:
        seg = get_segment_label(rec['read_id'])
        if seg is None:
            seg = "full"
        seg_groups.setdefault(seg, []).append(rec)

    best_segs = {}
    other_segs = {}
    one_diff_segs = {}
    tie_segs = {}

    for seg, recs in seg_groups.items():
        if len(recs) == 1:
            best_segs[seg] = recs
        else:
            # Sort by raw score to determine top two candidates
            sorted_recs_by_score = sorted(recs, key=similarity_score, reverse=True)

            # Check for one-difference condition
            if len(sorted_recs_by_score) >= 2:
                top1_rec = sorted_recs_by_score[0]
                top2_rec = sorted_recs_by_score[1]

                if abs(top1_rec['NM'] - top2_rec['NM']) == 1:
                    qcov1 = (top1_rec['qend'] - top1_rec['qstart']) / top1_rec['qlen'] if top1_rec['qlen'] > 0 else 0
                    qcov2 = (top2_rec['qend'] - top2_rec['qstart']) / top2_rec['qlen'] if top2_rec['qlen'] > 0 else 0
                    
                    if qcov1 >= 0.95 and qcov2 >= 0.95:
                        one_diff_segs[seg] = recs
                        continue
                    
                    cov1 = top1_rec['aln_block'] / top1_rec['tlen'] if top1_rec['tlen'] > 0 else 0
                    cov2 = top2_rec['aln_block'] / top2_rec['tlen'] if top2_rec['tlen'] > 0 else 0
                    rel_diff = abs(cov1 - cov2) / max(cov1, cov2) if max(cov1, cov2) > 0 else 0
                    if cov1 > 0 and cov2 > 0 and rel_diff < coverage_tolerance:
                        one_diff_segs[seg] = recs
                        continue

                # Check for inconclusive tie
                if similarity_score(top1_rec) == similarity_score(top2_rec):
                    tie_segs[seg] = recs
                    continue
            
            # If not one-diff or tie, select the single best and put others in 'other'
            best_rec_for_seg = select_best_alignment(recs) # This uses the standard select_best_alignment
            best_segs[seg] = [best_rec_for_seg]
            other_segs[seg] = [r for r in recs if r != best_rec_for_seg]
            
    return best_segs, other_segs, one_diff_segs, tie_segs

def process_best_matches(aln_list, coverage_tolerance=0.1):
    """
    Process the alignments for a read to get the best match from each segment group.
    Returns:
        - best_matches: list of best alignment records (one per segment group)
        - best_segs: dictionary keyed by segment label (for further use)
        - other_segs, one_diff_segs, tie_segs: as returned by process_read_segments.
    """
    best_segs, other_segs, one_diff_segs, tie_segs = process_read_segments(aln_list, coverage_tolerance)
    best_matches = []
    for seg, recs in best_segs.items():
        best_matches.extend(recs)
    return best_matches, best_segs, other_segs, one_diff_segs, tie_segs

def is_chimeric(groups, threshold=1000000):
    """
    Check for coordinate-based chimeric evidence between target groups (e.g., Col vs. Ler) for the full read.
    For each ecotype group, build a dictionary mapping segment labels to the best segment.
    For full-read alignments, use label 'full'.
    For each label common between groups, compare their target start positions.
    If any homologous pair shows an absolute difference greater than threshold, flag as chimeric.
    If no common label exists, fall back to comparing the best overall segment from each group.
    Returns (True, details) if chimeric; otherwise (False, details).
    """
    if "Col" not in groups or "Ler" not in groups:
        return (False, "Both Col and Ler groups not present")

    best_segments = {"Col": {}, "Ler": {}}
    for grp in ["Col", "Ler"]:
        for rec in groups[grp]:
            label = get_segment_label(rec['read_id'])
            if label is None:
                label = "full"
            aln_length = rec['qend'] - rec['qstart']
            
            current_best_for_label = best_segments[grp].get(label)
            if current_best_for_label is None or \
               similarity_score(rec) > similarity_score(current_best_for_label) or \
               (similarity_score(rec) == similarity_score(current_best_for_label) and aln_length > (current_best_for_label['qend'] - current_best_for_label['qstart'])):
                best_segments[grp][label] = rec

    comparisons = []
    chim_flag = False
    common_labels = set(best_segments["Col"].keys()).intersection(set(best_segments["Ler"].keys()))
    if common_labels:
        for label in common_labels:
            tstart_col = best_segments["Col"][label]['tstart']
            tstart_ler = best_segments["Ler"][label]['tstart']
            diff = abs(tstart_col - tstart_ler)
            comparisons.append(f"(Label {label}: Col tstart={tstart_col} vs Ler tstart={tstart_ler}, diff={diff} bp)")
            if diff > threshold:
                chim_flag = True
    else:
        if groups["Col"] and groups["Ler"]:
            # Need to pick the best overall alignment from each group to compare
            # This is done by taking the best alignment from the list of group alignments.
            # `select_best_alignment` here does not care about expected ecotype, just raw score.
            best_col = select_best_alignment(groups["Col"])
            best_ler = select_best_alignment(groups["Ler"])

            if best_col and best_ler:
                diff = abs(best_col['tstart'] - best_ler['tstart'])
                comparisons.append(f"(Overall: Col tstart={best_col['tstart']} vs Ler tstart={best_ler['tstart']}, diff={diff} bp)")
                if diff > threshold:
                    chim_flag = True
            else:
                details = "Could not find valid best alignments for comparison between Col and Ler groups."
                return (False, details)
        else:
            details = "One or both ecotype groups (Col/Ler) are empty for chimera comparison."
            return (False, details)

    details = "; ".join(comparisons)
    return (chim_flag, details)

def process_alignments(alignments, chimera_threshold=1000000, coverage_tolerance=0.1):
    """
    Process alignments (grouped by base read id) and classify each read into one exclusive category.
    Order of processing:
        1. Inversions (using all alignments)
        2. Determine the best match per segment group. If any best match is to the 'wrong' ecotype,
           classify the entire read as 'wrong_ecotype_best_alignment'.
        3. For remaining reads:
            - Check if any segment group qualifies as only_one_difference.
            - Check if any segment group is an inconclusive tie.
            - Check if the best matches map to different chromosomes (chimeric_different_chromosomes).
            - Else, check if they map to different arms (chimeric_different_arms).
            - Else, check for coordinate-based chimera.
        4. If none of the above, output the best matches as "best" (and others as "other").
    """
    best_output = {}
    other_output = {}
    one_diff_output = {}
    tie_output = {}
    inversions_output = {}
    chimeras_output = {}
    diff_arms_output = {}
    diff_chrom_output = {}
    wrong_ecotype_best_output = {} # New category

    for base_id, aln_list in alignments.items():
        print(f"Processing read {base_id} with {len(aln_list)} segments.")
        
        # 1. Inversions.
        if has_inversion(aln_list):
            inv_details = get_inversion_details(aln_list)
            print(f"Read {base_id} flagged as inversion. Details: {inv_details}")
            inversions_output[base_id] = [f"#{base_id}|inversion: {inv_details}\n"] + [rec['line'] for rec in aln_list]
            continue

        # 2. Get best matches per segment group based on raw score, and check for wrong ecotype mapping
        best_matches, best_segs, other_segs, one_diff_segs, tie_segs = process_best_matches(aln_list, coverage_tolerance)

        # Check for segments whose best alignment is to the "wrong" ecotype
        is_wrong_ecotype_mapping = False
        wrong_ecotype_details = []
        for seg_label, recs in best_segs.items():
            if recs: # Ensure there's a record in best_segs for this segment
                best_rec_for_seg = recs[0] # Get the single best record for this segment
                expected_ecotype = get_read_ecotype(best_rec_for_seg['read_id'])
                target_ecotype = "Col" if "Col" in best_rec_for_seg['target'] else ("Ler" if "Ler" in best_rec_for_seg['target'] else None)
                
                if expected_ecotype and target_ecotype and expected_ecotype != target_ecotype:
                    is_wrong_ecotype_mapping = True
                    wrong_ecotype_details.append(f"Segment {seg_label} (expected {expected_ecotype}) mapped best to {target_ecotype}")
                    
        if is_wrong_ecotype_mapping:
            print(f"Read {base_id} flagged as wrong_ecotype_best_alignment. Details: {'; '.join(wrong_ecotype_details)}")
            wrong_ecotype_best_output[base_id] = [f"#{base_id}|wrong_ecotype_best_alignment: {'; '.join(wrong_ecotype_details)}\n"] + \
                                                  [rec['line'] for rec in sorted(aln_list, key=lambda x: x['qstart'])]
            continue

        # 3a. Check one-difference.
        if one_diff_segs:
            lines = []
            for seg, recs in one_diff_segs.items():
                header = f"#{base_id}|only_one_difference:{seg}\n"
                sorted_recs = sorted(recs, key=lambda x: x['qstart'])
                lines.append(header)
                lines.extend([r['line'] for r in sorted_recs])
            one_diff_output[base_id] = lines
            continue

        # 3b. Check inconclusive tie.
        if tie_segs:
            lines = []
            for seg, recs in tie_segs.items():
                header = f"#{base_id}|inconclusive_tie:{seg}\n"
                sorted_recs = sorted(recs, key=lambda x: x['qstart'])
                lines.append(header)
                lines.extend([r['line'] for r in sorted_recs])
            tie_output[base_id] = lines
            continue

        # 3c. Check chimerism on best_matches.
        # First, check if best_matches map to different chromosomes.
        diff_chr_flag, diff_chr_details = is_different_chromosomes(best_matches)
        if diff_chr_flag:
            print(f"Read {base_id} flagged as chimeric_different_chromosomes. Details: {diff_chr_details}")
            sorted_best = sorted(best_matches, key=lambda x: x['qstart'])
            diff_chrom_output[base_id] = [f"#{base_id}|chimeric_different_chromosomes: {diff_chr_details}\n"] + [r['line'] for r in sorted_best]
            continue
        # Next, check if best_matches map to different arms.
        diff_arms_flag, diff_arms_details = is_different_arms_chimera(best_matches)
        if diff_arms_flag:
            print(f"Read {base_id} flagged as chimeric_different_arms. Details: {diff_arms_details}")
            sorted_best = sorted(best_matches, key=lambda x: x['qstart'])
            diff_arms_output[base_id] = [f"#{base_id}|chimeric_different_arms: {diff_arms_details}\n"] + [r['line'] for r in sorted_best]
            continue
        # Finally, check for coordinate-based chimera.
        groups = group_by_target(best_matches)
        chim_flag, chim_details = is_chimeric(groups, threshold=chimera_threshold)
        if chim_flag:
            print(f"Read {base_id} flagged as chimeric (coordinate-based). Details: {chim_details}")
            sorted_best = sorted(best_matches, key=lambda x: x['qstart'])
            chimeras_output[base_id] = [f"#{base_id}|chimeric: {chim_details}\n"] + [r['line'] for r in sorted_best]
            continue

        # 4. Otherwise, output best matches as "best".
        best_lines = []
        for seg, recs in best_segs.items():
            header = f"#{base_id}|best:{seg}\n"
            sorted_recs = sorted(recs, key=lambda x: x['qstart'])
            best_lines.append(header)
            best_lines.extend([r['line'] for r in sorted_recs])
        other_lines = []
        for seg, recs in other_segs.items():
            header = f"#{base_id}|other:{seg}\n"
            sorted_recs = sorted(recs, key=lambda x: x['qstart'])
            other_lines.append(header)
            other_lines.extend([r['line'] for r in sorted_recs])
        if best_lines:
            best_output[base_id] = best_lines
        if other_lines:
            other_output[base_id] = other_lines
        print(f"Read {base_id} classified; Best segments: {list(best_segs.keys())}; Other segments: {list(other_segs.keys())}")
    
    return best_output, other_output, one_diff_output, tie_output, inversions_output, chimeras_output, diff_arms_output, diff_chrom_output, wrong_ecotype_best_output

def write_output(filename, data_dict):
    """
    Write the dictionary of read blocks to a file.
    Each key's value is a list of lines.
    """
    with open(filename, "w") as f:
        for base_id in data_dict:
            for line in data_dict[base_id]:
                f.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Process PAF files to classify read segments based on mapping to Col vs. Ler."
    )
    parser.add_argument("paf_files", nargs="+", help="Input PAF files.")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory for result files.")
    parser.add_argument("-t", "--chimera_threshold", type=int, default=1000000,
                        help="Chimera threshold in bp (default 1,000,000 bp).")
    parser.add_argument("-c", "--coverage_tolerance", type=float, default=0.1,
                        help="Coverage tolerance for one-difference condition (default 0.1).")
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    test_read_id_parsing()
    
    alignments = load_alignments(args.paf_files)
    best_out, other_out, one_diff_out, tie_out, inversions_out, chimeras_out, diff_arms_out, diff_chrom_out, wrong_ecotype_best_out = process_alignments(
        alignments, chimera_threshold=args.chimera_threshold, coverage_tolerance=args.coverage_tolerance
    )

    write_output(os.path.join(outdir, "best_alignments.paf"), best_out)
    write_output(os.path.join(outdir, "other_alignments.paf"), other_out)
    write_output(os.path.join(outdir, "only_one_difference.paf"), one_diff_out)
    write_output(os.path.join(outdir, "inconclusive_tie.paf"), tie_out)
    write_output(os.path.join(outdir, "inversions.paf"), inversions_out)
    write_output(os.path.join(outdir, "chimeras.paf"), chimeras_out)
    write_output(os.path.join(outdir, "chimeric_different_arms.paf"), diff_arms_out)
    write_output(os.path.join(outdir, "chimeric_different_chromosomes.paf"), diff_chrom_out)
    write_output(os.path.join(outdir, "wrong_ecotype_best_alignment.paf"), wrong_ecotype_best_out) # New output file

    sys.stderr.write(f"Processed {len(alignments)} base read IDs.\n")
    sys.stderr.write(f"Best alignments written: {len(best_out)}\n")
    sys.stderr.write(f"Other alignments written: {len(other_out)}\n")
    sys.stderr.write(f"Only-one-difference reads written: {len(one_diff_out)}\n")
    sys.stderr.write(f"Inconclusive tie reads written: {len(tie_out)}\n")
    sys.stderr.write(f"Inversions written: {len(inversions_out)}\n")
    sys.stderr.write(f"Coordinate-based chimeras written: {len(chimeras_out)}\n")
    sys.stderr.write(f"Chimeric_different_arms reads written: {len(diff_arms_out)}\n")
    sys.stderr.write(f"Chimeric_different_chromosomes reads written: {len(diff_chrom_out)}\n")
    sys.stderr.write(f"Wrong ecotype best alignment reads written: {len(wrong_ecotype_best_out)}\n")

if __name__ == "__main__":
    main()
