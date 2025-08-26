#!/bin/bash
script="$1"
base_dir="$2"

declare -A patterns=(
    ["mm_lr_hqae_alnTo"]="lrhqae"
    ["mm_lr_hq_alnTo"]="lrhq"
    ["mm_ont_alnTo"]="mm_ont"
    ["mm_sr_alnTo"]="mm_sr"
    ["wm_mapont_alnTo"]="wm"   
)

echo "Starting analysis in directory: ${base_dir}"
echo "Available files:"
ls -la "${base_dir}/" 2>/dev/null || echo "Directory not accessible"

for pattern in "${!patterns[@]}"; do
    echo "=========================================="
    echo "Processing pattern: ${pattern} -> ${patterns[$pattern]}"
    
    # Use nullglob to handle empty matches gracefully
    shopt -s nullglob
    input_files=("${base_dir}/"*${pattern}*)
    shopt -u nullglob
    
    # Check if any files were actually found
    if [ ${#input_files[@]} -eq 0 ]; then
        echo "INFO: No files found for pattern *${pattern}* in ${base_dir}"
        echo "This is normal if this alignment type wasn't performed."
        
        # Create empty output directory for consistency
        output_dir="${base_dir}/summary_${patterns[$pattern]}"
        mkdir -p "$output_dir"
        echo "Created empty output directory: $output_dir"
        continue
    fi
    
    # Verify files actually exist (additional safety check)
    actual_files=()
    for file in "${input_files[@]}"; do
        if [ -f "$file" ]; then
            actual_files+=("$file")
        fi
    done
    
    if [ ${#actual_files[@]} -eq 0 ]; then
        echo "INFO: No valid files found for pattern *${pattern}*"
        output_dir="${base_dir}/summary_${patterns[$pattern]}"
        mkdir -p "$output_dir"
        continue
    fi
    
    # Create output directory
    output_dir="${base_dir}/summary_${patterns[$pattern]}"
    mkdir -p "$output_dir"
    
    echo "SUCCESS: Found ${#actual_files[@]} files for pattern: $pattern"
    echo "Files found:"
    for file in "${actual_files[@]}"; do
        echo "  - $(basename "$file")"
    done
    
    echo "Running Python script..."
    # Call Python script with actual files
    python "$script" "${actual_files[@]}" -o "$output_dir"
    
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Python script failed for pattern $pattern (exit code: $exit_code)"
        echo "This might be due to empty files or data issues."
    else
        echo "SUCCESS: Completed processing for pattern $pattern"
    fi
    echo "Output directory: $output_dir"
done

echo "=========================================="
echo "Bash script completed successfully."
echo "Summary directories created:"
ls -la "${base_dir}"/summary_* 2>/dev/null || echo "No summary directories found"
