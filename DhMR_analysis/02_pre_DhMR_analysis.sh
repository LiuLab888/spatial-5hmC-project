#!/bin/sh

sh merge_bismark_bed_files.sh sample_name cluster_name CpG

cd .../results/DMR_analysis/sample_name/raw_data/Cluster_name

less merge_sorted_filter_end_new_chr.bed | awk -F "\t" '{$7=$5+$6; print $1,$2,$7,$5}' OFS="\t" > merge_sorted_filter_end_new_chr_dss.bed

sed -i '1 i chr\tpos\tN\tX' merge_sorted_filter_end_new_chr_dss.bed
