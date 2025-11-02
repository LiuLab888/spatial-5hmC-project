#!/bin/sh
usage="sh <script.sh> <index> <Output_path1> <Output_path2> <bismark> <bismark_index> "

if [ $# -ne 5 ]; then

        echo "$usage"

        exit

fi


j=$1

Output_path1=$2

Output_path2=$3

bismark=$4

bismark_index=$5


trim_galore --cores 4  -a CTATCTCTTATACACATCT  --clip_R1 43  --basename ${j}_spatial_crick --output_dir ${Output_path2} ${Output_path1}/${j}_spatial_Read1.extract_crick.fq.gz  


cd ${Output_path2}

#part two

${bismark}/bismark  --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_trimmed.fq.gz 


#1

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_spatial_crick_unmapped_1 --output_dir ${Output_path2}  ${j}_spatial_crick_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_spatial_crick_unmapped_1_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_spatial_crick_unmapped_1_trimmed.fq.gz

${bismark}/bismark --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_unmapped_1_trimmed.fq.gz


#2

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_spatial_crick_unmapped_2 --output_dir ${Output_path2}  ${j}_spatial_crick_unmapped_1_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_spatial_crick_unmapped_2_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_spatial_crick_unmapped_2_trimmed.fq.gz

${bismark}/bismark --pbat  --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_unmapped_2_trimmed.fq.gz


#3

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_spatial_crick_unmapped_3 --output_dir ${Output_path2}  ${j}_spatial_crick_unmapped_2_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_spatial_crick_unmapped_3_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_spatial_crick_unmapped_3_trimmed.fq.gz

${bismark}/bismark --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_unmapped_3_trimmed.fq.gz


#4

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_spatial_crick_unmapped_4 --output_dir ${Output_path2}  ${j}_spatial_crick_unmapped_3_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_spatial_crick_unmapped_4_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_spatial_crick_unmapped_4_trimmed.fq.gz

${bismark}/bismark --pbat  --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_unmapped_4_trimmed.fq.gz


#5

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_spatial_crick_unmapped_5 --output_dir ${Output_path2}  ${j}_spatial_crick_unmapped_4_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_spatial_crick_unmapped_5_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_spatial_crick_unmapped_5_trimmed.fq.gz

${bismark}/bismark  --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_spatial_crick_unmapped_5_trimmed.fq.gz




