use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashSet;
use regex::Regex;
use rayon::prelude::*;

#[derive(Debug)]
struct ProfileEntry {
    code: String,
    quality: i32,
}

fn create_curated_profile(raw: &str, readmer_cutoff: i32) -> Vec<ProfileEntry> {
    raw.split(',')
        .filter_map(|entry| {
            let parts: Vec<&str> = entry.split(':').collect();
            if parts.len() == 2 {
                if let Ok(quality) = parts[1].parse::<i32>() {
                    if quality >= readmer_cutoff && parts[0] != "0" {
                        return Some(ProfileEntry {
                            code: parts[0].to_string(),
                            quality,
                        });
                    }
                }
            }
            None
        })
        .collect()
}

fn create_supercurated_profile(curated: &[ProfileEntry], threshold: i32) -> String {
    if curated.is_empty() {
        return String::new();
    }

    let mut supercurated_groups = Vec::new();
    let mut current_code = &curated[0].code;
    let mut count = 1;

    for entry in &curated[1..] {
        if entry.code == *current_code {
            count += 1;
        } else {
            if count >= threshold {
                supercurated_groups.push(format!("{}{}", count, current_code));
            }
            current_code = &entry.code;
            count = 1;
        }
    }

    if count >= threshold {
        supercurated_groups.push(format!("{}{}", count, current_code));
    }

    supercurated_groups.join(",")
}

fn create_ultracurated_profile(supercurated_str: &str) -> String {
    if supercurated_str.is_empty() {
        return String::new();
    }

    let re = Regex::new(r"(\d+)([A-Z0-9]+)").unwrap();
    let mut curated = Vec::new();

    for entry in supercurated_str.split(',') {
        if let Some(caps) = re.captures(entry) {
            let count: i32 = caps[1].parse().unwrap_or(0);
            let code = caps[2].to_string();
            curated.push((count, code));
        }
    }

    if curated.is_empty() {
        return String::new();
    }

    let mut ultra_curated_groups = Vec::new();
    let mut current_code = curated[0].1.clone();
    let mut total_count = curated[0].0;

    for (count, code) in &curated[1..] {
        if *code == current_code {
            total_count += count;
        } else {
            ultra_curated_groups.push(format!("{}{}", total_count, current_code));
            current_code = code.clone();
            total_count = *count;
        }
    }

    ultra_curated_groups.push(format!("{}{}", total_count, current_code));
    ultra_curated_groups.join(",")
}

fn hybrid_label(ultracurated_profile: &str) -> String {
    if ultracurated_profile.is_empty() {
        return "no|0".to_string();
    }

    let re = Regex::new(r"(\d+)([A-Z0-9]+)").unwrap();
    let mut distinct_codes = HashSet::new();

    for group in ultracurated_profile.split(',') {
        if let Some(caps) = re.captures(group) {
            distinct_codes.insert(caps[2].to_string());
        }
    }

    let count_codes = distinct_codes.len();
    if count_codes >= 2 {
        format!("yes|{}", count_codes)
    } else {
        format!("no|{}", count_codes)
    }
}

fn classify_hybrid(ultracurated_profile: &str) -> String {
    if ultracurated_profile.is_empty() {
        return "None".to_string();
    }

    let re = Regex::new(r"(\d+)([A-Z0-9]+)").unwrap();
    let mut parsed = Vec::new();

    for group in ultracurated_profile.split(',') {
        if let Some(caps) = re.captures(group) {
            let count: i32 = caps[1].parse().unwrap_or(0);
            let code = caps[2].to_string();
            parsed.push((count, code));
        }
    }

    if parsed.is_empty() {
        return "None".to_string();
    }

    let codes: Vec<&String> = parsed.iter().map(|(_, code)| code).collect();
    let distinct_codes: HashSet<&String> = codes.iter().cloned().collect();

    if distinct_codes.len() != 2 {
        return "Unknown".to_string();
    }

    let n = parsed.len();
    if n == 2 {
        return "SCO".to_string();
    }

    // Check alternating pattern
    let alternating = codes.windows(2).all(|w| w[0] != w[1]);
    if !alternating {
        return "Unknown".to_string();
    }

    if n == 3 && codes[0] == codes[2] {
        return "NCO".to_string();
    }

    if n == 4 {
        return "GC-CO".to_string();
    }

    "Unknown".to_string()
}

fn ectopic_classification(code1: &str, code2: &str) -> String {
    let re = Regex::new(r"([A-Z]{2})(\d+)").unwrap();
    
    let caps1 = re.captures(code1);
    let caps2 = re.captures(code2);
    
    if let (Some(c1), Some(c2)) = (caps1, caps2) {
        let prefix1 = &c1[1];
        let prefix2 = &c2[1];
        let num1: i32 = c1[2].parse().unwrap_or(0);
        let num2: i32 = c2[2].parse().unwrap_or(0);
        
        // Homologous check
        if &code1[1..] == &code2[1..] {
            return "no".to_string();
        }
        
        // Edge case
        if code1.chars().next() == code2.chars().next() && num1 == num2 && prefix1 != prefix2 {
            return "no|edge".to_string();
        }
        
        // Ectopic classification
        let access = if code1.chars().next() == code2.chars().next() {
            "intraaccession"
        } else {
            "interaccession"
        };
        
        let chrom = if num1 == num2 {
            "intrachromosome"
        } else {
            "interchromosome"
        };
        
        format!("yes|{}|{}", access, chrom)
    } else {
        "None".to_string()
    }
}

fn process_read(readname: &str, raw_profile: &str, thresholds: &[i32], readmer_cutoff: i32) -> Vec<(String, i32, String, String, String, String, String, String)> {
    let curated = create_curated_profile(raw_profile, readmer_cutoff);
    let curated_str = curated.iter()
        .map(|e| format!("{}:{}", e.code, e.quality))
        .collect::<Vec<_>>()
        .join(",");
    
    let mut results = Vec::new();
    
    for &threshold in thresholds {
        let supercurated = create_supercurated_profile(&curated, threshold);
        let ultra_curated = create_ultracurated_profile(&supercurated);
        let hybrid = hybrid_label(&ultra_curated);
        
        let (hybrid_type, ectopic) = if hybrid.starts_with("yes|2") {
            let ht = classify_hybrid(&ultra_curated);
            
            // Extract distinct codes for ectopic classification
            let re = Regex::new(r"\d+([A-Z0-9]+)").unwrap();
            let codes: HashSet<String> = ultra_curated.split(',')
                .filter_map(|group| re.captures(group))
                .map(|caps| caps[1].to_string())
                .collect();
            
            let ectopic = if codes.len() == 2 {
                let codes_vec: Vec<String> = codes.into_iter().collect();
                ectopic_classification(&codes_vec[0], &codes_vec[1])
            } else {
                "None".to_string()
            };
            
            (ht, ectopic)
        } else {
            ("None".to_string(), "None".to_string())
        };
        
        results.push((
            readname.to_string(),
            threshold,
            curated_str.clone(),
            supercurated,
            ultra_curated,
            hybrid,
            hybrid_type,
            ectopic,
        ));
    }
    
    results
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 || args.len() > 5 {
        eprintln!("Usage: {} <input_file> <output_file> <readmer_cutoff> [num_threads]", args[0]);
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];
    let readmer_cutoff: i32 = args[3].parse().expect("Invalid readmer_cutoff");

    // Set number of threads (default to all available cores)
    if args.len() == 5 {
        let num_threads: usize = args[4].parse().expect("Invalid num_threads");
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .unwrap();
        eprintln!("Using {} threads for parallel processing", num_threads);
    }

    let thresholds = [10, 20, 30, 35, 40, 45, 50, 75, 100];

    // CHUNK-BASED PROCESSING: Process in chunks to limit memory usage
    // Chunk size: number of reads to process at once (adjust based on available memory)
    const CHUNK_SIZE: usize = 50000; // Process 50k reads at a time (~500MB-1GB per chunk)

    eprintln!("Processing input file in chunks of {} reads...", CHUNK_SIZE);

    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let output = File::create(output_file)?;
    let mut writer = BufWriter::new(output);

    // Write header
    writeln!(writer, "readname\tthreshold\tcurated_profile\tsuper_curated_profile\tultra_curated_profile\thybrid\thybrid_type\tectopic")?;

    let mut total_processed = 0usize;
    let mut chunk_number = 0usize;

    loop {
        // Read a chunk of reads
        let mut read_chunk = Vec::with_capacity(CHUNK_SIZE);

        for _ in 0..CHUNK_SIZE {
            if let Some(line) = lines.next() {
                let line = line?;

                if line.starts_with('>') {
                    let readname = line[1..].trim().to_string();

                    if let Some(profile_line) = lines.next() {
                        let raw_profile = profile_line?;
                        read_chunk.push((readname, raw_profile));
                    }
                }
            } else {
                break; // End of file
            }
        }

        if read_chunk.is_empty() {
            break; // No more reads to process
        }

        chunk_number += 1;
        let chunk_size = read_chunk.len();
        eprintln!("Processing chunk {} ({} reads)...", chunk_number, chunk_size);

        // Process this chunk in parallel
        let chunk_results: Vec<Vec<(String, i32, String, String, String, String, String, String)>> =
            read_chunk.par_iter()
                .map(|(readname, raw_profile)| {
                    process_read(readname, raw_profile, &thresholds, readmer_cutoff)
                })
                .collect();

        // Write chunk results immediately (stream to disk, don't accumulate)
        for read_results in chunk_results {
            for row in read_results {
                writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        row.0, row.1, row.2, row.3, row.4, row.5, row.6, row.7)?;
            }
        }

        total_processed += chunk_size;
        eprintln!("Completed {} reads total", total_processed);

        // Clear the chunk to free memory before next iteration
        drop(read_chunk);
    }

    writer.flush()?;
    eprintln!("TSV file '{}' created successfully.", output_file);
    eprintln!("Total reads processed: {}", total_processed);

    Ok(())
}
