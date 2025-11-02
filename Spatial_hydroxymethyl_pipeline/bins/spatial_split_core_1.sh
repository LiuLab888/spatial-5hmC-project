#!/bin/sh
usage="sh <script.sh> <name> <index> <work_path>"

if [ $# -ne 3 ]; then 
	echo "$usage"
	exit
fi

name=$1
index=$2
work_path=$3


less  ${work_path}/${name}_spatial_2.extract.barcode.order | awk -F "\t" -v barcode=${index} '{if($1==barcode) print $2}' > ${index}.txt

less  ${work_path}/${name}_spatial_1.extract_index.txt  | awk 'NR==FNR{a[FNR]=$0}NR!=FNR{print a[$1]}' - ${index}.txt  >  ${index}_fastq.txt

less ${index}_fastq.txt | awk -F " " '{$3=substr($1,2,length($1)-2);print $0}' OFS=" " | awk -F " " '{$4=$3$2; print $4}' > ${index}_fastq_subtract_1.txt

rm -fr ${index}_fastq.txt ${index}.txt
