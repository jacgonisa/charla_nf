use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
struct Segment {
    read_name: String,
    start: i32,
    end: i32,
    region: String,
    token_count: i32,
}

fn process_kmer_profile(read_name: &str, profile: &str, k: i32) -> Vec<Segment> {
    let tokens: Vec<&str> = profile.split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    
    if tokens.is_empty() {
        return Vec::new();
    }
    
    // Extract only the region label from each token (ignoring any number after the colon)
    let labels: Vec<String> = tokens.iter()
        .map(|token| token.split(':').next().unwrap_or("").to_string())
        .collect();
    
    let mut segments = Vec::new();
    let mut current_label = labels[0].clone();
    let mut start_index = 0;
    
    for (i, label) in labels.iter().enumerate() {
        if label != &current_label {
            let block_start = start_index;
            let block_end = (i as i32 - 1) + k; // token (i-1) covers bases [i-1, i-1+k)
            let token_count = i as i32 - start_index;
            segments.push(Segment {
                read_name: read_name.to_string(),
                start: block_start,
                end: block_end,
                region: current_label.clone(),
                token_count,
            });
            current_label = label.clone();
            start_index = i as i32;
        }
    }
    
    // Append the final block
    let block_start = start_index;
    let block_end = (labels.len() as i32 - 1) + k;
    let token_count = labels.len() as i32 - start_index;
    segments.push(Segment {
        read_name: read_name.to_string(),
        start: block_start,
        end: block_end,
        region: current_label,
        token_count,
    });
    
    segments
}

fn merge_consecutive_by_code(segments: Vec<Segment>) -> Vec<Segment> {
    if segments.is_empty() {
        return Vec::new();
    }
    
    // Filter out isolated small segments
    let filtered: Vec<Segment> = segments.into_iter()
        .filter(|seg| !((seg.end - seg.start) < 50 && seg.token_count < 50))
        .collect();
    
    if filtered.is_empty() {
        return Vec::new();
    }
    
    let mut merged = Vec::new();
    let mut current = filtered[0].clone();
    
    for seg in &filtered[1..] {
        // If the new segment has the same region code, merge it with the current chain
        if seg.region == current.region {
            // Merge: keep the original start and update the end to the new segment's end
            current.end = seg.end;
            current.token_count += seg.token_count;
        } else {
            merged.push(current.clone());
            current = seg.clone();
        }
    }
    merged.push(current);
    merged
}

fn read_curated_tsv(tsv_file: &str) -> std::io::Result<HashMap<String, HashSet<String>>> {
    let file = File::open(tsv_file)?;
    let reader = BufReader::new(file);
    let mut curated = HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.is_empty() {
            continue;
        }
        
        let read_name = parts[0].to_string();
        if let Some(last_col) = parts.last() {
            let allowed_codes: HashSet<String> = last_col.split('-')
                .map(|s| s.to_string())
                .collect();
            curated.insert(read_name, allowed_codes);
        }
    }
    
    Ok(curated)
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <input_kmer_profile.fa> <curated_tsv_file> <output.bed> <k_size>", args[0]);
        std::process::exit(1);
    }
    
    let kmer_profile_file = &args[1];
    let curated_tsv_file = &args[2];
    let outfile = &args[3];
    let k: i32 = args[4].parse().expect("Invalid k_size");
    
    // Read curated TSV file
    let curated = read_curated_tsv(curated_tsv_file)?;
    eprintln!("Loaded {} curated reads", curated.len());
    
    // Open input and output files
    let input_file = File::open(kmer_profile_file)?;
    let reader = BufReader::new(input_file);
    
    let output_file = File::create(outfile)?;
    let mut writer = BufWriter::new(output_file);
    
    let mut current_read: Option<String> = None;
    let mut profile_lines = Vec::new();
    let mut processed_count = 0;
    
    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with('>') {
            // Process the previous read, if any
            if let Some(read_name) = &current_read {
                if !profile_lines.is_empty() {
                    let profile = profile_lines.join("");
                    if let Some(allowed) = curated.get(read_name) {
                        let segments = process_kmer_profile(read_name, &profile, k);
                        // Keep only segments with allowed codes
                        let filtered_segments: Vec<Segment> = segments.into_iter()
                            .filter(|seg| allowed.contains(&seg.region))
                            .collect();
                        
                        let merged = merge_consecutive_by_code(filtered_segments);
                        for seg in merged {
                            writeln!(writer, "{}\t{}\t{}\t{}", 
                                    seg.read_name, seg.start, seg.end, seg.region)?;
                        }
                    }
                    
                    processed_count += 1;
                    if processed_count % 1000 == 0 {
                        eprintln!("Processed {} reads...", processed_count);
                    }
                }
            }
            
            current_read = Some(line[1..].trim().to_string());
            profile_lines.clear();
        } else {
            profile_lines.push(line.trim().to_string());
        }
    }
    
    // Process the final read
    if let Some(read_name) = &current_read {
        if !profile_lines.is_empty() {
            let profile = profile_lines.join("");
            if let Some(allowed) = curated.get(read_name) {
                let segments = process_kmer_profile(read_name, &profile, k);
                let filtered_segments: Vec<Segment> = segments.into_iter()
                    .filter(|seg| allowed.contains(&seg.region))
                    .collect();
                
                let merged = merge_consecutive_by_code(filtered_segments);
                for seg in merged {
                    writeln!(writer, "{}\t{}\t{}\t{}", 
                            seg.read_name, seg.start, seg.end, seg.region)?;
                }
            }
            processed_count += 1;
        }
    }
    
    writer.flush()?;
    eprintln!("Processed {} total reads", processed_count);
    println!("BED file '{}' created.", outfile);
    
    Ok(())
}
