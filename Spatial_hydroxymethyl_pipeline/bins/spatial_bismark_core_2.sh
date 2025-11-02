#!/bin/sh
usage="sh <script.sh> <index> <Output_path1> <Output_path2> <bismark> <bismark_index>"

if [ $# -ne 5 ]; then
 
	echo "$usage"

	exit

fi


j=$1

Output_path1=$2

Output_path2=$3

bismark=$4

bismark_index=$5


$bismark/deduplicate_bismark -s --bam ${Output_path1}/${j}.bam 

$bismark/bismark_methylation_extractor -s   --multicore 1  --bedGraph   --genome_folder  $bismark_index  ${j}.deduplicated.bam


rm -fr CpG*${j}*  CHH*${j}*  CHG*${j}*

rm -fr ${j}*.M-bias.txt  ${j}*_splitting_report.txt


samtools view ${j}.deduplicated.bam | less | awk -v barcode=${j} 'END{print barcode"\t"NR}' >> depth_matrix.txt


less ${j}.deduplicated.bismark.cov.gz | awk -F "\t"  '{if($1=="1" || $1=="2" || $1=="3" || $1=="4" || $1=="5" || $1=="6" || $1=="7" || $1=="8" || $1=="9" || $1=="10" || $1=="11" || $1=="12" || $1=="13" || $1=="14" || $1=="15" || $1=="16" || $1=="17" || $1=="18" || $1=="19" || $1=="X" || $1=="Y" || $1=="MT") print $0}' OFS="\t" | awk -F "\t" -v barcode=${j} 'BEGIN{C=0;mC=0}{if($5>0 || $6>0){mC=mC+$5;C=C+$5+$6}}END{print barcode"\t"mC/C}' >> methylation_matrix.txt


less ${j}.deduplicated.bismark.cov.gz | awk -F "\t"  '{if($1=="1" || $1=="2" || $1=="3" || $1=="4" || $1=="5" || $1=="6" || $1=="7" || $1=="8" || $1=="9" || $1=="10" || $1=="11" || $1=="12" || $1=="13" || $1=="14" || $1=="15" || $1=="16" || $1=="17" || $1=="18" || $1=="19" || $1=="X" || $1=="Y" || $1=="MT") print $0}' OFS="\t" | awk -F "\t" -v barcode=${j} 'BEGIN{C=0}{if($5>0 || $6>0){C=C+1}}END{print barcode"\t"C}' >> coverage_matrix.txt


