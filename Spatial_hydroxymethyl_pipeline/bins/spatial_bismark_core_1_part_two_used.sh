#!/bin/sh
usage="sh <script.sh> <bismark> <cpu> <index> <output_prefix> <fqgz>"

if [ $# -ne 5 ]; then

        echo "$usage"

        exit

fi

bismark=$1
bismark_index=$2
index=$bismark_index
fqgz=$3
j=$4
prefix=$j
cpu=$5

trim_galore --cores 4  -a CTATCTCTTATACACATCT  --clip_R1 43  --basename `basename ${prefix}` $fqgz


#part two

${bismark}/bismark  --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_trimmed.fq.gz


#1

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_1  ${j}_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_unmapped_1_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_1_trimmed.fq.gz

${bismark}/bismark --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_unmapped_1_trimmed.fq.gz


#2

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_2  ${j}_unmapped_1_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_unmapped_2_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_2_trimmed.fq.gz

${bismark}/bismark --pbat  --unmapped --dovetail -p 8 $bismark_index  ${j}_unmapped_2_trimmed.fq.gz


#3

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_3  ${j}_unmapped_2_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_unmapped_3_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_3_trimmed.fq.gz

${bismark}/bismark --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_unmapped_3_trimmed.fq.gz


#4

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_4 ${j}_unmapped_3_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_unmapped_4_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_4_trimmed.fq.gz

${bismark}/bismark --pbat  --unmapped --dovetail -p 8 $bismark_index  ${j}_unmapped_4_trimmed.fq.gz


#5

trim_galore --cores 4  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_5 ${j}_unmapped_4_trimmed.fq.gz_unmapped_reads.fq.gz

mv ${j}_unmapped_5_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_5_trimmed.fq.gz

${bismark}/bismark  --pbat --unmapped --dovetail -p 8 $bismark_index  ${j}_unmapped_5_trimmed.fq.gz




