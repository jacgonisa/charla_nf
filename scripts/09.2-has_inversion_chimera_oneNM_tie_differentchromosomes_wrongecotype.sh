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

for pattern in "${!patterns[@]}"; do
    # Expand wildcard pattern to actual files
    input_files=("${base_dir}/"*_${pattern}_*)
    
    # Skip if no files found
    if [ ${#input_files[@]} -eq 0 ]; then
        echo "Warning: No files found for pattern: *_${pattern}_*"
        continue
    fi

    # Create output directory
    output_dir="${base_dir}/summary_${patterns[$pattern]}"
    mkdir -p "$output_dir"

    echo "Processing ${#input_files[@]} files for pattern: $pattern"
    python "$script" "${input_files[@]}" -o "$output_dir"
done
