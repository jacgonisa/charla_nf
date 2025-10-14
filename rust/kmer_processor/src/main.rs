use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::sync::Arc;

use clap::{Arg, Command};
use glob::glob;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;

// Custom error type for better error handling
#[derive(Debug)]
struct ProcessError(String);

impl std::fmt::Display for ProcessError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ProcessError {}

impl From<std::io::Error> for ProcessError {
    fn from(err: std::io::Error) -> Self {
        ProcessError(format!("IO error: {}", err))
    }
}

impl From<glob::GlobError> for ProcessError {
    fn from(err: glob::GlobError) -> Self {
        ProcessError(format!("Glob error: {}", err))
    }
}

impl From<glob::PatternError> for ProcessError {
    fn from(err: glob::PatternError) -> Self {
        ProcessError(format!("Pattern error: {}", err))
    }
}

impl From<&str> for ProcessError {
    fn from(err: &str) -> Self {
        ProcessError(err.to_string())
    }
}

impl From<String> for ProcessError {
    fn from(err: String) -> Self {
        ProcessError(err)
    }
}

#[derive(Debug, Clone)]
struct ReadCount {
    read_name: String,
    count: usize,
}

#[derive(Debug)]
struct FileResult {
    filename: String,
    reads: HashMap<String, usize>,
}

fn main() -> Result<(), ProcessError> {
    let matches = Command::new("K-mer Count Processor")
        .version("1.0")
        .about("Process k-mer count files into a count matrix")
        .arg(
            Arg::new("base_dir")
                .long("base-dir")
                .value_name("DIR")
                .help("Directory containing *.txt files")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .value_name("FILE")
                .help("Output TSV file path")
                .required(true),
        )
        .arg(
            Arg::new("total_reads")
                .long("total-reads")
                .value_name("NUM")
                .help("Total number of reads per file")
                .required(true)
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(
            Arg::new("cpus")
                .long("cpus")
                .value_name("NUM")
                .help("Number of CPU cores to use")
                .default_value("1")
                .value_parser(clap::value_parser!(usize)),
        )
        .get_matches();

    let base_dir = matches.get_one::<String>("base_dir").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let total_reads_per_file = *matches.get_one::<usize>("total_reads").unwrap();
    let cpus = *matches.get_one::<usize>("cpus").unwrap();

    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(output_file).parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Find all .txt files (both directly in base_dir and in subdirectories)
    let pattern_direct = format!("{}/*.txt", base_dir);
    let pattern_recursive = format!("{}/**/*.txt", base_dir);

    let mut files: Vec<PathBuf> = Vec::new();

    // Try direct pattern first
    if let Ok(entries) = glob(&pattern_direct) {
        files.extend(entries.filter_map(Result::ok));
    }

    // Then try recursive pattern
    if let Ok(entries) = glob(&pattern_recursive) {
        files.extend(entries.filter_map(Result::ok));
    }

    // Remove duplicates (in case a file matches both patterns)
    files.sort();
    files.dedup();

    if files.is_empty() {
        return Err(ProcessError(format!("No *.txt files found in {}", base_dir)));
    }

    println!("[INFO] Found {} *.txt files in: {}", files.len(), base_dir);

    // Sort files using custom logic
    files.sort_by(|a, b| custom_sort_key(a).cmp(&custom_sort_key(b)));
    
    let filenames: Vec<String> = files
        .iter()
        .map(|f| f.file_name().unwrap().to_string_lossy().to_string())
        .collect();
    
    println!("[INFO] Sorted files: {:?}", filenames);

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cpus)
        .build_global()
        .unwrap();

    println!("[INFO] Using {} threads for parallel processing", cpus);

    // Process files in parallel with enhanced progress tracking
    let start_time = Instant::now();
    let multi_progress = Arc::new(MultiProgress::new());
    
    // Create master progress bar for overall file processing
    let master_pb = multi_progress.add(ProgressBar::new(files.len() as u64));
    master_pb.set_style(
        ProgressStyle::default_bar()
            .template("🚀 Overall Progress: [{elapsed_precise}] {bar:40.green/red} {pos}/{len} files completed ({percent}%) ETA: {eta}")
            .map_err(|e| ProcessError(format!("Master progress bar template error: {}", e)))?
            .progress_chars("█▉▊▋▌▍▎▏  ")
    );
    master_pb.set_message("Processing k-mer files...");
    
    // Clone master_pb for use in the parallel processing
    let master_pb_clone = master_pb.clone();
    
    let results: Result<Vec<FileResult>, ProcessError> = files
        .par_iter()
        .enumerate()
        .map(|(index, file_path)| {
            let result = process_file(
                file_path, 
                total_reads_per_file, 
                Arc::clone(&multi_progress),
                index + 1,  // 1-indexed for display
                files.len()
            );
            
            // Update master progress bar when file completes
            if result.is_ok() {
                master_pb_clone.inc(1);
            }
            
            result
        })
        .collect();
    
    master_pb.finish_with_message("✅ All files processed successfully!");
    let results = results?;

    println!("[INFO] Processing completed in {:.2}s", start_time.elapsed().as_secs_f64());
    println!("[📊] Summary:");
    println!("    • Files processed: {}/20", results.len());
    println!("    • Average time per file: {:.2}s", start_time.elapsed().as_secs_f64() / results.len() as f64);
    println!("    • Total throughput: {:.0} reads/second", 
             (results.len() * total_reads_per_file) as f64 / start_time.elapsed().as_secs_f64());

    // Combine results
    let mut all_reads: HashMap<String, HashMap<String, usize>> = HashMap::new();
    
    for file_result in results {
        println!("[✓] Processed: {}", file_result.filename);
        for (read_name, count) in file_result.reads {
            all_reads
                .entry(read_name)
                .or_insert_with(HashMap::new)
                .insert(file_result.filename.clone(), count);
        }
    }

    // Write output
    write_output_table(&all_reads, &filenames, output_file)?;
    
    println!("[✔] Output written to: {}", output_file);
    println!("[INFO] Processed {} unique reads across {} files", 
             all_reads.len(), files.len());

    Ok(())
}

fn custom_sort_key(path: &Path) -> (usize, usize, String) {
    let filename = path.file_name().unwrap().to_string_lossy();
    
    let acc_order = if filename.contains("Col") { 1 }
                   else if filename.contains("Ler") { 2 }
                   else { 0 };
    
    let reg_order = if filename.contains("CEN") { 4 }
                   else if filename.contains("ARMS") { 3 }
                   else { 0 };
    
    let chr_num: String = filename.chars().filter(|c| c.is_ascii_digit()).collect();
    
    (acc_order, reg_order, chr_num)
}

fn process_file(
    file_path: &Path,
    total_reads: usize,
    multi_progress: Arc<MultiProgress>,
    file_number: usize,
    total_files: usize,
) -> Result<FileResult, ProcessError> {
    let filename = file_path.file_name().unwrap().to_string_lossy().to_string();
    
    // Create progress bar for this specific file
    let pb = multi_progress.add(ProgressBar::new(total_reads as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&format!("📁 File {}/{}: [{{elapsed_precise}}] {{bar:30.cyan/blue}} {{pos:>8}}/{{len:8}} {{msg}}", file_number, total_files))
            .map_err(|e| ProcessError(format!("Progress bar template error: {}", e)))?
            .progress_chars("█▉▊▋▌▍▎▏  "),
    );
    pb.set_message(format!("{}", filename));

    let file = File::open(file_path)?;
    let reader = BufReader::with_capacity(128 * 1024, file); // Increased buffer to 128KB
    
    let mut reads: HashMap<String, usize> = HashMap::with_capacity(total_reads); // Pre-allocate capacity
    let mut lines = reader.lines();
    let batch_size = 15000; // Slightly larger batches
    let mut batch = Vec::with_capacity(batch_size);
    let mut processed_count = 0;

    loop {
        // Read batch of header-count pairs
        batch.clear();
        
        for _ in 0..batch_size {
            // Read header line
            if let Some(header_line) = lines.next() {
                let header = header_line?;
                if !header.starts_with('>') {
                    continue;
                }
                
                // Read count line
                if let Some(count_line) = lines.next() {
                    let counts = count_line?;
                    batch.push((header, counts));
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        if batch.is_empty() {
            break;
        }

        // Process batch in parallel
        let batch_results: Vec<ReadCount> = batch
            .par_iter()
            .filter_map(|(header, count_line)| {
                parse_read_count(header, count_line).ok()
            })
            .collect();

        // Update results and progress
        for result in batch_results {
            reads.insert(result.read_name, result.count);
        }
        
        processed_count += batch.len();
        pb.set_position(processed_count.min(total_reads) as u64);
        
        // Optional: Update message with current read count
        if processed_count % 50000 == 0 {  // Update every 50k reads
            pb.set_message(format!("{} ({:.1}%)", filename, (processed_count as f64 / total_reads as f64) * 100.0));
        }
    }

    pb.finish_with_message(format!("✅ {} ({} reads)", filename, reads.len()));

    Ok(FileResult { filename, reads })
}

fn parse_read_count(header: &str, count_line: &str) -> Result<ReadCount, ProcessError> {
    // Parse read name from header
    let read_name = header
        .trim()
        .split_whitespace()
        .next()
        .ok_or_else(|| ProcessError("Invalid header format".to_string()))?
        .trim_start_matches('>')
        .to_string();

    // Parse counts and calculate non-zero count
    let count_line = count_line.trim().trim_end_matches(',');
    let non_zero_count = count_line
        .split(',')
        .filter_map(|s| s.parse::<i32>().ok())
        .filter(|&x| x != 0)
        .count();

    Ok(ReadCount {
        read_name,
        count: non_zero_count,
    })
}

fn write_output_table(
    all_reads: &HashMap<String, HashMap<String, usize>>,
    filenames: &[String],
    output_file: &str,
) -> Result<(), ProcessError> {
    let mut output = File::create(output_file)?;
    
    // Write header
    write!(output, "ReadName")?;
    for filename in filenames {
        write!(output, "\t{}", filename)?;
    }
    writeln!(output)?;
    
    // Sort read names for consistent output
    let mut read_names: Vec<_> = all_reads.keys().collect();
    read_names.sort();
    
    // Write data rows
    for read_name in read_names {
        write!(output, "{}", read_name)?;
        
        let read_data = &all_reads[read_name];
        for filename in filenames {
            let count = read_data.get(filename).unwrap_or(&0);
            write!(output, "\t{}", count)?;
        }
        writeln!(output)?;
    }
    
    Ok(())
}
