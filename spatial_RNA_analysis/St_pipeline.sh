#!/bin/bash

# Validate input arguments
if [ $# -lt 11 ]; then
    echo "Error: Missing required arguments!"
    echo "Usage: $0 <threads> <memory> <sample_name> <output_path> <read1> <read2> <STAR_index> <annotation_gtf> <contaminant_index> <barcode_file> <gene_model>"
    exit 1
fi

# Assign input arguments
ppn=$1                # Number of threads (e.g., 12)
mem=$2                # Memory in GB (e.g., 60)
name=$3               # Sample name (e.g., E18_0910)
path=$4               # Output directory (e.g., /path/to/output)
FW=$5                 # Forward read file (e.g., /path/to/sample_2.extract.fq.gz)
RV=$6                 # Reverse read file (e.g., /path/to/sample_1.extract.fq.gz)
MAP=$7                # STAR index (e.g., /path/to/STAR_index)
ANN=$8                # Annotation GTF file (e.g., /path/to/annotation.gtf)
CONT=$9               # Contaminant index (e.g., /path/to/rDNA_STAR_index)
ID_File=${10}         # Barcode reference file (e.g., /path/to/combine_barcode.txt)
gene_model=${11}      # Gene model file for RSeQC (e.g., /path/to/gene_model.bed)

# Define output directories
OUTPUT="$path/1_stpipeline"
TMP="$OUTPUT/tmp"

# Create necessary directories
mkdir -p "$TMP"

# Run Spatial Transcriptomics Pipeline
st_pipeline_run.py \
  --output-folder "$OUTPUT" \
  --temp-folder "$TMP" \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --ids "$ID_File" \
  --ref-map "$MAP" \
  --ref-annotation "$ANN" \
  --expName "$name" \
  --htseq-no-ambiguous \
  --verbose \
  --threads "$ppn" \
  --log-file "$OUTPUT/${name}_log.txt" \
  --star-two-pass-mode \
  --no-clean-up \
  --contaminant-index "$CONT" \
  --disable-clipping \
  --min-length-qual-trimming 15 \
  --star-sort-mem-limit "${mem}000000000"  \
  "$FW" "$RV"

# Convert Ensembl IDs to Gene Symbols
convertEnsemblToNames.py \
  --annotation "$ANN" \
  --output "$OUTPUT/${name}_stdata.symbol.tsv" \
  "$OUTPUT/${name}_stdata.tsv"

# Quality Assessment (Optional)
# Uncomment to enable barcode quality assessment
# st_qa-new.py --input-data "$OUTPUT/${name}_stdata.tsv"
# st_qa.py --input-data "$OUTPUT/${name}_stdata.tsv"

# Generate BAM Index for Visualization
samtools index "$TMP/mapped.bam"

# Generate BigWig Coverage File (Optional)
# Uncomment the following line to enable BAM coverage calculation
# bamCoverage -b "$TMP/mapped.bam" --normalizeUsing RPKM -o "$TMP/mapped.fpkm.bw" -p "$ppn" --binSize 10

# Run Read Distribution Analysis (Optional)
# Uncomment the following line to check read distribution
# read_distribution.py -i "$TMP/mapped.bam" -r "$gene_model"

# Additional Read Counts (Optional)
# Uncomment to count mapped reads
# echo -e "\n## Reads after removing rRNA ##"
# samtools view -c "$TMP/contaminated_clean.bam"
# echo -e "## Reads after mapping ##"
# samtools view -c "$TMP/mapped.bam"

echo "Pipeline completed successfully!"
