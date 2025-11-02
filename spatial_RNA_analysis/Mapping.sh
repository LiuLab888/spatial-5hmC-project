#!/bin/sh

# Define core parameters
core=$core
mem=$mem
name=$name
filepath=$filepath  # Base directory for input/output

# Define input/output file paths
RNA1=${filepath}/${name}_1.fq.gz
RNA2=${filepath}/${name}_2.fq.gz
RNA_R1=${filepath}/${name}_1.extract.fq.gz
RNA_R2=${filepath}/${name}_2.extract.fq.gz
out_path=${filepath}

# Define reference paths (Replace these with generalized paths)
STAR_INDEX_DIR=/path/to/STAR_index
GTF_FILE=/path/to/genome_annotation.gtf
RDNA_INDEX=/path/to/rDNA_STAR_index
BARCODE_FILE=/path/to/barcode_list.txt
GENE_MODEL=/path/to/gene_model.bed

# Create the SLURM job script
cat >> 00.spR_analysis.${name}.sh << 'EOF'
#!/bin/sh
#SBATCH --job-name=spR_${name}      # Job name
#SBATCH --cpus-per-task=$core       # Number of CPU cores
#SBATCH --mem=${mem}G               # Memory allocation
#SBATCH --output=${out_path}/${name}.log  # Log file

# Load necessary environments
export PATH="/home/user/anaconda3/bin:$PATH"
source /home/user/.bashrc
conda activate Env.UMItools

# Change to output directory
cd $out_path

# Step 1: Extract barcodes
/path/to/Extract_barcode.sh $name $RNA1 $RNA2 $out_path/

# Extract barcodes from read 2
zcat $RNA_R2 | awk 'NR%4==2{print substr($1,1,16)}' > $out_path/${name}_2.extract.barcode

# Step 2: Check barcode quality
Rscript /path/to/Check_barcode.R -b $BARCODE_FILE -p $out_path -n ${name}_2.extract.barcode

# Cleanup intermediate barcode file
rm $out_path/${name}_2.extract.barcode

# Step 3: Run spatial transcriptomics pipeline
conda deactivate
conda activate my_st_env

/path/to/st_pipeline.sh $core $mem $name $out_path \
$RNA_R2 $RNA_R1 $STAR_INDEX_DIR $GTF_FILE $RDNA_INDEX $BARCODE_FILE $GENE_MODEL
EOF

# Submit the job script
sbatch 0.spR_analysis.${name}.sh
