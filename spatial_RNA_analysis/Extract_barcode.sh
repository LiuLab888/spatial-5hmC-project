#!/bin/bash

# Check if at least four arguments are provided
if [ $# -lt 4 ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <sample_name> <RNA_R1> <RNA_R2> <output_path>"
    exit 1
fi

# Assign input arguments to variables
name=$1       # Sample name (e.g., E18.5_stain_20um_RNA)
RNA_R1=$2     # Read 1 input FASTQ file
RNA_R2=$3     # Read 2 input FASTQ file
out_path=$4   # Output directory

# Create output directory if it does not exist
mkdir -p ${out_path}

# Step 1: Extract UMIs from Read 1
umi_tools extract --extract-method=regex \
--bc-pattern2="^(?P<umi_1>.)(?P<discard_1>ATCGGCGTACGACT){s<=1}.{8}(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}.{8}(?P<discard_3>CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG){s<=2,i<=1,d<=1}.{10}(?P<discard_4>.*)$" \
-I ${RNA_R1} \
-S ${out_path}/${name}_1.extract.fq.gz \
--read2-in=${RNA_R2} \
--read2-out=${out_path}/${name}_2.extract.fq.gz \
-L ${out_path}/extract.log

# Print completion message
echo "UMI extraction completed. Output saved in ${out_path}."
