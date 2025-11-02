#!/bin/sh
usage="sh <script.sh> <name> <index> <work_path> <split>"

if [ $# -ne 4 ]; then 
	echo "$usage"
	exit
fi

name=$1
index=$2
work_path=$3
k=$4


touch ${index}.sam

less ${work_path}/${name}_bismark_index${k} | awk 'NR==FNR{a[$1]=$2}NR!=FNR{print a[$1]}' -  ${index}_fastq_subtract_1.txt | grep -v '^$' > ${index}_bismark_index.order${k}

less ${work_path}/${name}.SE.merge.filtered${k} | awk 'NR==FNR{a[NR]=$0}NR!=FNR{print a[$1]}' - ${index}_bismark_index.order${k} >> ${index}.sam
