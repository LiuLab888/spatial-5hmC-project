#!/bin/sh
usage="sh <script.sh> <name> <cluster> <base>"

if [ $# -ne 3 ];then

        echo "$usage"

        exit

fi


name=$1

cluster=$2

base=$3


combine_barcode=.../merge_step/${name}_new/Cluster_${cluster}_cell_index_one_row.txt

work_path=.../results/A_preWork_outputdir_${name}_merge_${base}

mkdir -p ${work_path}/bismark/matrix/all_Clusters/Cluster${cluster}

Output_path=${work_path}/bismark/matrix/all_Clusters/Cluster${cluster}


cd ${work_path}




touch ${Output_path}/merge.bed


for  j  in `cat ${combine_barcode}`

do

zcat  bismark/matrix/part_*/\${j}.deduplicated.bismark.cov.gz  >> ${Output_path}/merge.bed

done


sort -k1,1 -k2,2n ${Output_path}/merge.bed > ${Output_path}/merge_sorted.bed

less ${Output_path}/merge_sorted.bed |  awk -F "\t"  '{if(\$1=="1" || \$1=="2" || \$1=="3" || \$1=="4" || \$1=="5" || \$1=="6" || \$1=="7" || \$1=="8" || \$1=="9" || \$1=="10" || \$1=="11" || \$1=="12" || \$1=="13" || \$1=="14" || \$1=="15" || \$1=="16" || \$1=="17" || \$1=="18" || \$1=="19" || \$1=="X" || \$1=="Y" || \$1=="MT") print \$0}' OFS="\t" > ${Output_path}/merge_sorted_filter.bed



less ${Output_path}/merge_sorted_filter.bed | awk -F "\t" '{ print \$1":"\$2":"\$3"\t"\$4"\t"\$5"\t"\$6}' | awk -F "\t" '{a[\$1] += \$2; b[\$1] += \$3; c[\$1] += \$4; num[\$1]++}END {for(i in a){print(i, a[i]/num[i], b[i], c[i])}}' OFS="\t" | sort -k1,1n > ${Output_path}/merge_sorted_filter_end_new.bed


sed -i 's/:/\t/g' ${Output_path}/merge_sorted_filter_end_new.bed


less  ${Output_path}/merge_sorted_filter_end_new.bed | awk -F "\t" '{\$1="chr"\$1; print \$0}' OFS="\t"  > ${Output_path}/merge_sorted_filter_end_new_chr.bed

