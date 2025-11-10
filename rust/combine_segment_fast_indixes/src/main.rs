use crate::io::BufWriter;
use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, Seek, SeekFrom, Write, ErrorKind},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    time::Instant,
};
use csv::ReaderBuilder;

/// CLI arguments
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Root directory containing hapmer count files (e.g., <parent>_<region>_<chr>_k<kmer>.txt)
    #[clap(long)]
    input_dir: PathBuf,

    /// Tab-delimited k-mer presence/absence matrix (TSV)
    #[clap(long)]
    kmer_tsv: PathBuf,

    /// QC file: FASTA-style headers with comma-separated quality scores per read
    #[clap(long)]
    qc_file: PathBuf,

    /// Where to write the final hybrid profiles
    #[clap(long)]
    output_file: PathBuf,

    /// Directory to store or load JSON indexes
    #[clap(long)]
    index_dir: PathBuf,

    /// Force rebuilding all indexes even if they exist
    #[clap(long, action)]
    rebuild_index: bool,

    /// Minimum k-mer count to include a read
    #[clap(long, default_value_t = 50)]
    min_count: u32,

    /// Verbose logging
    #[clap(short, long, action)]
    verbose: bool,
}

/// In-memory representation of the k-mer TSV
#[derive(Debug, Serialize, Deserialize)]
struct KmerIndex {
    reads: Vec<String>,
    headers: Vec<String>,
    matrix: Vec<Vec<u32>>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Verbose logging macro
    macro_rules! vlog {
        ($($arg:tt)+) => {
            if args.verbose {
                eprintln!($($arg)+)
            }
        }
    }

    let total_start = Instant::now();

    // Ensure the index directory exists
    fs::create_dir_all(&args.index_dir)
        .with_context(|| format!("creating {}", args.index_dir.display()))?;

    // Paths to each JSON index
    let kmer_json = args.index_dir.join("kmer_index.json");
    let qc_json = args.index_dir.join("qc_index.json"); // Will now store offsets
    let cam_json = args.index_dir.join("cam_index.json"); // offsets-only
    let abbrev_json = args.index_dir.join("abbrev_map.json");

    // Determine if indexes need rebuilding
    let needs_rebuild = args.rebuild_index
        || !kmer_json.exists()
        || !qc_json.exists()
        || !cam_json.exists()
        || !abbrev_json.exists();

    let kmer_index: KmerIndex;
    let qc_index: HashMap<String, u64>; // Changed to store offset
    let abbrev_map: HashMap<String, String>;
    let cambridge_filenames: Vec<String>;
    let cambridge_files: Vec<PathBuf>;
    let cam_index: HashMap<String, HashMap<String, u64>>;

    if needs_rebuild {
        vlog!("🔧 Rebuilding JSON indexes in {}…", args.index_dir.display());

        // 1) Build k-mer index
        let start = Instant::now();
        vlog!("    Loading k-mer presence matrix from {}…", args.kmer_tsv.display());
        let (reads, headers, matrix) = load_kmer_presence(&args.kmer_tsv)?;
        kmer_index = KmerIndex { reads: reads.clone(), headers: headers.clone(), matrix };
        serde_json::to_writer_pretty(File::create(&kmer_json)?, &kmer_index)?;
        vlog!("    ✅ Wrote {} ({} reads, {} files) in {:?}",
            kmer_json.display(), kmer_index.reads.len(), kmer_index.headers.len() - 1, start.elapsed());

        // Prepare list of reads from k-mer index for filtering other indexes
        let kmer_reads_set: HashSet<String> = kmer_index.reads.iter().cloned().collect();
        vlog!("    Relevant reads identified from k-mer matrix: {}", kmer_reads_set.len());


        // 2) Build QC offset index
        let start = Instant::now();
        vlog!("    Building QC offset index for relevant reads from {}…", args.qc_file.display());
        qc_index = build_qc_offset_index(&args.qc_file, &kmer_reads_set, args.verbose)?;
        serde_json::to_writer_pretty(File::create(&qc_json)?, &qc_index)?;
        vlog!("    ✅ Wrote {} ({} reads indexed) in {:?}", qc_json.display(), qc_index.len(), start.elapsed());

        // 3) Build abbreviation map
        let start = Instant::now();
        vlog!("    Building abbreviation map…");
        abbrev_map = kmer_index.headers.iter().skip(1).filter_map(|fname| {
            make_abbrev(fname).map(|ab| (fname.clone(), ab))
        }).collect();
        serde_json::to_writer_pretty(File::create(&abbrev_json)?, &abbrev_map)?;
        vlog!("    ✅ Wrote {} ({} files) in {:?}", abbrev_json.display(), abbrev_map.len(), start.elapsed());

        // 4) Build hapmer file offset index
        cambridge_filenames = kmer_index.headers.iter().skip(1).cloned().collect();

        // NEW (Corrected)
        cambridge_files = cambridge_filenames.iter()
        .map(|fname| {
        // Simply join the input_dir directly with the filename
        args.input_dir.join(fname) // <-- Corrected: No subdirectories assumed
        })
        .collect();


        let start = Instant::now();
        vlog!("    Building hapmer file offset index for relevant reads…");
        cam_index = load_or_build_cam_index(&cambridge_files, &kmer_reads_set, &cam_json, args.verbose)?;
        vlog!("    ✅ Wrote {} (offsets for {} files) in {:?}",
            cam_json.display(), cam_index.len(), start.elapsed());

    } else {
        vlog!("📂 JSON indexes found—skipping rebuild.\n");

        // ─────────────────────────────────────────────────────────
        // Load JSON indexes back into memory if not rebuilt
        vlog!("📥 Loading k-mer index from {}…", kmer_json.display());
        let start = Instant::now();
        kmer_index = serde_json::from_reader(File::open(&kmer_json)?)?;
        vlog!("    ✅ Loaded k-mer index ({} reads, {} files) in {:?}",
            kmer_index.reads.len(), kmer_index.headers.len() - 1, start.elapsed());

        vlog!("📥 Loading QC offset index from {}…", qc_json.display());
        let start = Instant::now();
        qc_index = serde_json::from_reader(File::open(&qc_json)?)?;
        vlog!("    ✅ Loaded QC offset index ({} reads) in {:?}", qc_index.len(), start.elapsed());

        vlog!("📥 Loading abbreviation map from {}…", abbrev_json.display());
        let start = Instant::now();
        abbrev_map = serde_json::from_reader(File::open(&abbrev_json)?)?;
        vlog!("    ✅ Loaded abbreviation map ({} files) in {:?}", abbrev_map.len(), start.elapsed());

        vlog!("📥 Loading hapmer file offset index from {}…", cam_json.display());
        let start = Instant::now();
        cam_index = serde_json::from_reader(File::open(&cam_json)?)?;
        vlog!("    ✅ Loaded hapmer file offset index (offsets for {} files) in {:?}", cam_index.len(), start.elapsed());

        // Need filenames/paths even if loading indexes
        cambridge_filenames = kmer_index.headers.iter().skip(1).cloned().collect();
        cambridge_files = cambridge_filenames.iter()
            .map(|fname| {
                let parts: Vec<&str> = fname.split('_').collect();
                let ec = parts.get(0).unwrap_or(&"");
                let region = if parts.get(1).map_or(false, |r| r.contains("CEN")) { "CEN" } else { "ARMS" };
                args.input_dir.join(ec).join(region).join(fname)
            })
            .collect();
        // ─────────────────────────────────────────────────────────
    }

    // Filter reads by min_count
    vlog!("🔍 Filtering reads based on min_count ≥ {}…", args.min_count);
    let start = Instant::now();
    let filtered_reads: Vec<String> = kmer_index
        .reads
        .par_iter()
        .zip(&kmer_index.matrix)
        .filter_map(|(r, row)| row.iter().any(|&c| c >= args.min_count).then(|| r.clone()))
        .collect();
    vlog!("📈 {} reads pass min_count ≥ {} (filtered from {} in {:?})",
        filtered_reads.len(), args.min_count, kmer_index.reads.len(), start.elapsed());

    // Map read → row index
    let idx_map: HashMap<String, usize> = kmer_index
        .reads
        .into_iter()
        .enumerate()
        .map(|(i, r)| (r, i))
        .collect();

    // Prepare output writer
    let out_f = OpenOptions::new()
        .create(true).truncate(true).write(true)
        .open(&args.output_file)?;
    let writer = Arc::new(Mutex::new(BufWriter::new(out_f)));

    // Progress bar for processing reads
    let pb = ProgressBar::new(filtered_reads.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("Processing Reads [{bar:40}] {pos}/{len} ({eta})")?
            .progress_chars("=>-"),
    );

    // ─────────────────────────────────────────────────────────
    // Parallel processing loop
    vlog!("🧬 Generating hybrid profiles…");
    let processing_start = Instant::now();
    filtered_reads
        .par_iter()
        .progress_with(pb.clone())
        .for_each(|read| {
            if let Some(&row_i) = idx_map.get(read) {
                let row = &kmer_index.matrix[row_i];
                // Initialize hybrid from QC: "0:score"
                let mut hybrid: Vec<String> = if let Some(&qc_offset) = qc_index.get(read) {
                    match read_qc_scores_at(&args.qc_file, qc_offset) {
                        Ok(scores) => scores.iter().map(|q| format!("0:{}", q)).collect(),
                        Err(e) => {
                            if args.verbose {
                                eprintln!("Warning: Error reading QC scores for read {} at offset {}: {}. Using NA scores.", read, qc_offset, e);
                            }
                            vec!["NA".into(); cambridge_files.len()] // Fallback length
                        }
                    }
                } else {
                    if args.verbose {
                        eprintln!("Warning: Read {} found in k-mer index but not QC index. Using NA scores.", read);
                    }
                    vec!["NA".into(); cambridge_files.len()] // Fallback length
                };

                // For each hapmer file column
                for (col_i, &count) in row.iter().enumerate() {
                    if count < args.min_count {
                        continue;
                    }
                    let fname = &cambridge_filenames[col_i];
                    if let Some(file_map) = cam_index.get(fname) {
                        if let Some(&offset) = file_map.get(read) {
                            // Pass read_id and verbose flag for debugging in read_cam_profile_at
                            match read_cam_profile_at(&cambridge_files[col_i], offset, read, args.verbose) {
                                Ok(profile) => {
                                    if profile.len() == hybrid.len() {
                                        if let Some(ab) = abbrev_map.get(fname) {
                                            for i in 0..profile.len() {
                                                if profile[i] != "0" && hybrid[i].starts_with("0:") {
                                                    let qc = hybrid[i].splitn(2, ':').nth(1).unwrap_or("NA");
                                                    hybrid[i] = format!("{}:{}", ab, qc);
                                                }
                                            }
                                        }
                                    } else {
                                        // Warning for mismatch if profile was read but wrong length
                                        if args.verbose {
                                            eprintln!("Warning: Profile length mismatch for read {} in file {}. Expected {}, got {}. Skipping profile merge for this file.",
                                                read, fname, hybrid.len(), profile.len());
                                        }
                                    }
                                }
                                Err(e) => {
                                    if args.verbose {
                                        eprintln!("Warning: Error reading profile for read {} at offset {} in file {}: {}. Skipping profile merge for this file.",
                                            read, offset, cambridge_files[col_i].display(), e);
                                    }
                                }
                            }
                        }
                    }
                }

                // Write out the hybrid profile
                let line = format!(">{}\n{}\n", read, hybrid.join(","));
                if let Ok(mut w) = writer.lock() {
                    let _ = w.write_all(line.as_bytes());
                } else {
                    if args.verbose {
                        eprintln!("Warning: Failed to acquire writer lock for read {}", read);
                    }
                }
            } else {
                if args.verbose {
                    eprintln!("Warning: Read {} not found in k-mer index map. Skipping.", read);
                }
            }
        });
    pb.finish_with_message("🎉 Done processing reads");
    vlog!("⏱️ Read processing completed in {:?}", processing_start.elapsed());

    vlog!("Overall execution time: {:?}", total_start.elapsed());

    Ok(())
}

/// Build or load the **offset** index for all hapmer files.
/// Only indexes reads specified in the candidates HashSet.
/// Writes/reads a JSON: HashMap< filename, HashMap< read, byte_offset > >.
/// Includes a progress bar for file processing.
/// **Refined to explicitly find the profile line after a header.**
fn load_or_build_cam_index(
    files: &[PathBuf],
    candidates: &HashSet<String>,
    idx_file: &Path,
    verbose: bool,
) -> io::Result<HashMap<String, HashMap<String, u64>>> {

    if idx_file.exists() {
        if verbose {
            eprintln!("⏳ Loading hapmer file index from {} …", idx_file.display());
        }
        let start = Instant::now();
        let f = File::open(idx_file)?;
        let m: HashMap<String, HashMap<String, u64>> = serde_json::from_reader(f)?;
        if verbose {
            eprintln!("    ✅ Loaded offsets for {} files in {:?}", m.len(), start.elapsed());
        }
        return Ok(m);
    }

    if verbose {
        eprintln!("🔨 Building hapmer file offset index …");
    }

    let start = Instant::now();
    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("Indexing hapmer files [{bar:40}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let map = files
        .par_iter()
        .progress_with(pb)
        .filter_map(|file| {
            match File::open(file) {
                Ok(f) => {
                    let mut reader = BufReader::new(f);
                    let mut line = String::new();
                    let mut m = HashMap::new();

                    loop {
                        // Get position before reading the current line
                        let current_line_start_pos = match reader.stream_position() {
                            Ok(pos) => pos,
                            Err(_) => {
                                if verbose { eprintln!("    ❌ Error getting stream position in {}", file.display()); }
                                return None; // Skip this file
                            }
                        };

                        match reader.read_line(&mut line) {
                            Ok(bytes) => {
                                if bytes == 0 { // EOF
                                    break;
                                }
                                // line contains the data, reader is positioned after it.
                            }
                            Err(_) => { // Error reading line
                                if verbose { eprintln!("    ❌ Error reading line in {}", file.display()); }
                                return None; // Skip this file on error
                            }
                        }

                        if line.starts_with('>') {
                            let key = line[1..].trim().to_string();
                            if candidates.contains(&key) {

                                // Found a relevant header. Now find the *next* non-empty, non-header line.
                                // This assumes the profile line is the first such line after the header.
                                let mut found_profile_line = false;
                                loop {
                                    let pos_before_reading_next = match reader.stream_position() {
                                        Ok(pos) => pos,
                                        Err(_) => {
                                            if verbose { eprintln!("    ❌ Error getting stream position while searching for profile after '{}' in {}", key, file.display()); }
                                            return None; // Skip this file
                                        }
                                    };
                                    let mut temp_line = String::new();
                                    match reader.read_line(&mut temp_line) {
                                        Ok(bytes) => {
                                            if bytes == 0 { // EOF while searching for profile
                                                if verbose {
                                                    eprintln!("    ⚠️  EOF while searching for profile after header '{}' in file {}", key, file.display());
                                                }
                                                // Did not find profile line before EOF
                                                break;
                                            }
                                            // Check if this line is the profile line
                                            let trimmed_temp_line = temp_line.trim();
                                            if !trimmed_temp_line.is_empty() && !trimmed_temp_line.starts_with('>') {
                                                // Found what we assume is the profile line
                                                m.insert(key.clone(), pos_before_reading_next); // Store its position
                                                found_profile_line = true;
                                                break; // Found profile, stop inner loop
                                            }
                                            // else: it was an empty line or another header, continue searching
                                        }
                                        Err(_) => {
                                            if verbose { eprintln!("    ❌ Error reading line while searching for profile after '{}' in file {}", key, file.display()); }
                                            return None; // Skip this file on error
                                        }
                                    }
                                } // End of inner loop searching for profile line

                                // After the inner loop, reader's position is after the line that terminated the loop
                                // (either the found profile line, or the last line before EOF/error).

                                if !found_profile_line && verbose {
                                    // This case is already handled by the inner loop's EOF break or error,
                                    // but this makes it explicit that we didn't index it.
                                    // eprintln!("    ⚠️  No valid profile line found after header '{}' in file {}", key, file.display());
                                }
                            }
                        }
                        line.clear(); // Clear line for next read cycle in outer loop
                    } // End of outer loop reading lines
                    Some((file.file_name()?.to_string_lossy().into_owned(), m))
                },
                Err(e) => {
                    if verbose {
                        eprintln!("    ❌ Error opening file {}: {}", file.display(), e);
                    }
                    None
                }
            }
        })
        .collect::<HashMap<_, _>>();

    let f = File::create(idx_file)?;
    serde_json::to_writer_pretty(f, &map)?;
    if verbose {
        eprintln!("    ✅ Saved hapmer file index to {} in {:?}", idx_file.display(), start.elapsed());
    }
    Ok(map)
}

/// Build a HashMap from relevant read headers to the byte offset of their score line in the QC file.
/// Only stores offsets for reads present in the relevant_reads HashSet.
/// Verbose mode will print diagnostics to stderr.
fn build_qc_offset_index(
    path: &Path,
    relevant_reads: &HashSet<String>,
    verbose: bool,
) -> io::Result<HashMap<String, u64>> {
    if verbose {
        eprintln!("🔍 Building QC offset index from {}…", path.display());
    }

    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let mut qc_offset_map = HashMap::new();
    let mut line = String::new();

    let total = relevant_reads.len() as u64;
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("Indexing QC file    [{bar:40}] {pos}/{len} reads ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    loop {
        let header_start = reader.stream_position()?;
        match reader.read_line(&mut line) {
            Ok(0) => break,              // EOF
            Ok(_) => {},
            Err(e) => {
                eprintln!("⚠ Warning: error reading QC file {}: {}. Stopping.", path.display(), e);
                break;
            }
        }

        if line.starts_with('>') {
            let read_id = line[1..].trim().to_string();
            line.clear();

            if relevant_reads.contains(&read_id) {
                // Capture the offset at the very start of the **scores** line:
                let score_offset = reader.stream_position()?;

                match reader.read_line(&mut line) {
                    Ok(bytes) if bytes > 0 => {
                        qc_offset_map.insert(read_id.clone(), score_offset);
                        pb.inc(1);
                        if verbose {
                            eprintln!(
                                "  • Indexed read {}: header at {}, scores at {}; first scores bytes: {}",
                                read_id,
                                header_start,
                                score_offset,
                                &line.trim()[..line.trim().find(',').unwrap_or(10).min(line.trim().len())]
                            );
                        }
                    }
                    Ok(_) => {
                        eprintln!(
                            "⚠ Warning: EOF immediately after header '{}' in {}",
                            read_id,
                            path.display()
                        );
                    }
                    Err(e) => {
                        eprintln!(
                            "⚠ Warning: error reading scores for '{}' in {}: {}",
                            read_id,
                            path.display(),
                            e
                        );
                    }
                }
            } else {
                // skip scores line for irrelevant read
                let mut dummy = String::new();
                let _ = reader.read_line(&mut dummy);
            }
        }
        line.clear();
    }

    pb.finish_with_message("✅ QC index built");
    if verbose {
        eprintln!(
            "🔖 Finished QC indexing: {} of {} reads indexed",
            qc_offset_map.len(),
            total
        );
    }

    Ok(qc_offset_map)
}

/// Jump to offset in file (QC file) and read the line containing comma-separated scores.
fn read_qc_scores_at(file: &Path, offset: u64) -> io::Result<Vec<String>> {
    let mut f = File::open(file)?;
    f.seek(SeekFrom::Start(offset))?;
    let mut reader = BufReader::new(f);
    let mut scores_line = String::new();
    reader.read_line(&mut scores_line)?;
    let scores: Vec<String> = scores_line.trim().split(',').map(String::from).collect();
    Ok(scores)
}

/// Load the k-mer presence/absence matrix from a TSV file.
/// Returns a tuple containing:
/// - Vec<String>: Read identifiers (first column).
/// - Vec<String>: Header row (including the first column label).
/// - Vec<Vec<u32>>: The binary presence/absence matrix (as u32).
fn load_kmer_presence(path: &Path) -> Result<(Vec<String>, Vec<String>, Vec<Vec<u32>>)> {
    let f = File::open(path)
        .with_context(|| format!("failed to open k-mer TSV: {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut csv_reader = ReaderBuilder::new().delimiter(b'\t').from_reader(reader);

    let headers = csv_reader.headers()?
        .iter()
        .map(String::from)
        .collect::<Vec<String>>();

    let mut reads = Vec::new();
    let mut matrix = Vec::new();

    for record in csv_reader.records() {
        let record = record?;
        if record.is_empty() {
            continue;
        }
        let read_id = record.get(0)
            .ok_or_else(|| anyhow::anyhow!("Missing read ID in k-mer TSV"))?
            .to_string();
        reads.push(read_id);

        let row = record.iter().skip(1)
            .map(|s| s.parse::<u32>().unwrap_or(0)) // Default to 0 on parse error
            .collect::<Vec<u32>>();
        matrix.push(row);
    }

    Ok((reads, headers, matrix))
}

/// Makes a shorter abbreviation from a hapmer filename.
/// Examples:
///   - Col_CEN_Chr1_k21.txt -> CC1
///   - C57BL6J_CHR_chr1_k21.txt -> CH1
/// Creates a three-character abbreviation from hapmer filename:
///   <parent>_<region>_<chr>_…  →  <P><R><chrom>
/// where P = first letter of parent name (uppercase),
///       R = "C" for CEN, "A" for ARMS, or "H" for CHR (centromere-unaware mode),
/// and  chrom = the number after "Chr" or "chr".
fn make_abbrev(filename: &str) -> Option<String> {
    let parts: Vec<&str> = filename.split('_').collect();
    let ecotype = parts.get(0)?;           // e.g. "Col", "Ler", "C57BL6J", "DBA2J"
    let region  = parts.get(1)?;           // e.g. "ARMS", "CEN", or "CHR"
    let chr_part = parts.get(2)?;          // e.g. "Chr1", "chr1", etc.

    // 1) First letter of the ecotype, uppercased
    let e = ecotype.chars().next()?.to_ascii_uppercase();

    // 2) Region: C for CEN, A for ARMS, H for CHR
    let r = if region.contains("CEN") {
        'C'
    } else if region.contains("ARMS") {
        'A'
    } else if region.contains("CHR") {
        'H'
    } else {
        // fallback or skip if unexpected
        return None;
    };

    // 3) Strip the "Chr" or "chr" prefix (case-insensitive) to get the chromosome number
    let chrom = chr_part
        .trim_start_matches("Chr")
        .trim_start_matches("chr");

    Some(format!("{}{}{}", e, r, chrom))
}





/// Reads a single hapmer profile at the given byte offset.
/// Expects the line at the offset to be the comma-separated profile.
/// Returns a Vec<String> of the profile values.
/// Handles potential errors during file access or reading.
fn read_cam_profile_at(file_path: &Path, offset: u64, read_id: &str, verbose: bool) -> Result<Vec<String>> {
    let file = File::open(file_path)
        .with_context(|| format!("Error opening hapmer file: {}", file_path.display()))?;
    let mut reader = BufReader::new(file);
    reader.seek(SeekFrom::Start(offset))
        .with_context(|| format!("Error seeking to offset {} in {}", offset, file_path.display()))?;

    let mut profile_line = String::new();
    reader.read_line(&mut profile_line)
        .with_context(|| format!("Error reading profile line at offset {} in {}", offset, file_path.display()))?;

    let profile: Vec<String> = profile_line.trim().split(',').map(String::from).collect();
    Ok(profile)
}


