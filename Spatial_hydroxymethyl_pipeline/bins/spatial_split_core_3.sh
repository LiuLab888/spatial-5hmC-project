#!/bin/sh
usage="sh <script.sh> <name> <index> <work_path>"

if [ $# -ne 3 ]; then 
	echo "$usage"
	exit
fi

name=$1
index=$2
work_path=$3


less  ${index}.sam | grep -v '^$' >  ${index}_final.sam

cat ${work_path}/header.txt ${index}_final.sam > ${index}_header.sam 

samtools view -bS ${index}_header.sam > ${index}.bam

