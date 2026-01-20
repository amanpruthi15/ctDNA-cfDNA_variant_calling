#!/bin/bash

# This script processes a list of genomic positions against BAM files in the current directory.
# It uses `samtools mpileup` and `awk` to extract base and mapping quality metrics for each position.
# The output is a tab-delimited table with columns for Sample, Position, Depth, Base counts (A, C, G, T), 
# Average Base Quality, and Average Mapping Quality.

# Usage:
#   ./calculate_base_proportion.sh <positions_file>
# 
# Arguments:
#   <positions_file>: A file containing genomic positions in the format "chr:start-end" (one per line).
#                     For example:
#                     chr1:1000-1001
#                     chr2:2000-2001
# 
# Example:
#   ./calculate_base_proportion.sh positions.txt

# Check if the input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <positions_file>"
    exit 1
fi

# File containing the positions
positions_file="$1"

# Print the table header
echo -e "Sample\tPosition\tDepth\tA\tC\tG\tT\tAvg-Base-Quality\tAvg-Mapping-Quality"

# Loop through each BAM file in the current directory
for bam_file in *.bam; do
    # Check if there are no BAM files
    if [ ! -e "$bam_file" ]; then
        echo "No BAM files found in the current directory."
        exit 1
    fi
    
    # Loop through each line in the positions file
    while IFS=$'\t' read -r position; do
        # Extract position components
        chr=$(echo "$position" | cut -d':' -f1)
        pos=$(echo "$position" | cut -d':' -f2 | cut -d'-' -f1)
        
        # Use samtools mpileup and awk to process the data
        samtools mpileup -r "$position" "$bam_file" 2>/dev/null | awk -v bam_file="$bam_file" -v position="$position" '
        {
            if ($1 ~ /^#/) next;  # Skip comment lines
            depth = $4;
            bases = tolower($5);
            qualities = $6;

            A = gsub("a", "a", bases);
            C = gsub("c", "c", bases);
            G = gsub("g", "g", bases);
            T = gsub("t", "t", bases);

            sum_base_quality = 0;

            for (i = 1; i <= length(qualities); i++) {
                sum_base_quality += (ord(substr(qualities, i, 1)) - 33);
            }

            avg_base_quality = (depth > 0) ? sum_base_quality / depth : "N/A";

            # Collect mapping qualities using samtools view
            cmd = "samtools view -q 0 " bam_file " " $1 ":" $2 "-" $2;
            total_mapping_quality = 0;
            mapping_quality_count = 0;

            while ((cmd | getline line) > 0) {
                split(line, fields, "\t");
                mapping_quality = fields[5];
                total_mapping_quality += mapping_quality;
                mapping_quality_count++;
            }
            close(cmd);

            avg_mapping_quality = (mapping_quality_count > 0) ? total_mapping_quality / mapping_quality_count : "N/A";

            printf "%s\t%s\t%d\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%s\t%s\n", bam_file, position, depth, (depth > 0 ? (A / depth) * 100 : 0), (depth > 0 ? (C / depth) * 100 : 0), (depth > 0 ? (G / depth) * 100 : 0), (depth > 0 ? (T / depth) * 100 : 0), avg_base_quality, avg_mapping_quality;
        }

        function ord(c) {
            return sprintf("%d", c);
        }'
    done < "$positions_file"
done
