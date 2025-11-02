#!/bin/bash
usage="sh <script.sh> <name> <Input_path> <Input_path1> "

if [ $# -ne 3 ]; then

    echo "$usage"

    exit
fi


name=$1

Input_path=$2

Input_path1=$3

#Read2 -- crick strand

cd ${Input_path}

zcat ${name}_spatial_2.extract_all_crick.fq.gz | less | awk 'NR%4==1{print $0}' > ${name}_spatial_2.extract_all_crick.filter.txt

less ${name}_spatial_2.extract_all_crick.filter.txt | cut -f1 -d " " | cut -f1 -d "_" > ${name}_spatial_2.extract_all_crick.filter_1_col.txt

rm -fr  ${name}_spatial_2.extract_all_crick.filter.txt



zcat ${Input_path1}/${name}_R2.fq.gz | less | awk 'NR%4==1{print $0}' >  ${name}_spatial_2.extract_all.filter1.txt

less ${name}_spatial_2.extract_all.filter1.txt  | cut -f1 -d " " > ${name}_spatial_2.extract_all.filter1_1_col.txt 

rm -fr ${name}_spatial_2.extract_all.filter1.txt



less ${name}_spatial_2.extract_all.filter1_1_col.txt | awk 'NR==FNR{a[$1]=NR}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_crick.filter_1_col.txt | grep -v '^$' > ${name}_spatial_2.extract_all_crick.filter.order1.txt

rm -fr  ${name}_spatial_2.extract_all.filter1_1_col.txt


less ${name}_spatial_2.extract_all_crick.filter.order1.txt | awk -F "\t" '{$2=($1-1)*4+1; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order1_1.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order1.txt | awk -F "\t" '{$2=($1-1)*4+2; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order1_2.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order1.txt | awk -F "\t" '{$2=($1-1)*4+3; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order1_3.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order1.txt | awk -F "\t" '{$2=($1-1)*4+4; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order1_4.txt


cat ${name}_spatial_2.extract_all_crick.filter.order1_1.txt  ${name}_spatial_2.extract_all_crick.filter.order1_2.txt  ${name}_spatial_2.extract_all_crick.filter.order1_3.txt  ${name}_spatial_2.extract_all_crick.filter.order1_4.txt | sort -k1,1n > ${name}_spatial_2.extract_all_crick.filter.order1_all.txt

rm -fr ${name}_spatial_2.extract_all_crick.filter.order1_1.txt ${name}_spatial_2.extract_all_crick.filter.order1_2.txt ${name}_spatial_2.extract_all_crick.filter.order1_3.txt ${name}_spatial_2.extract_all_crick.filter.order1_4.txt

zcat  ${Input_path1}/${name}_R2.fq.gz | less | awk 'NR==FNR{a[NR]=$0}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_crick.filter.order1_all.txt | gzip  > ${name}_spatial_Read2.extract_crick.fq.gz

rm -fr ${name}_spatial_2.extract_all_crick.filter.order1_all.txt



#Read1 -- crick strand


zcat ${Input_path1}/${name}_R1.fq.gz | less | awk 'NR%4==1{print $0}' >  ${name}_spatial_1.extract_all.filter2.txt

less ${name}_spatial_1.extract_all.filter2.txt  | cut -f1 -d " "  > ${name}_spatial_1.extract_all.filter2_1_col.txt 

rm -fr ${name}_spatial_1.extract_all.filter2.txt



less ${name}_spatial_1.extract_all.filter2_1_col.txt | awk 'NR==FNR{a[$1]=NR}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_crick.filter_1_col.txt | grep -v '^$' > ${name}_spatial_2.extract_all_crick.filter.order2.txt

rm -fr ${name}_spatial_1.extract_all.filter2_1_col.txt ${name}_spatial_2.extract_all_crick.filter_1_col.txt



less ${name}_spatial_2.extract_all_crick.filter.order2.txt | awk -F "\t" '{$2=($1-1)*4+1; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order2_1.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order2.txt | awk -F "\t" '{$2=($1-1)*4+2; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order2_2.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order2.txt | awk -F "\t" '{$2=($1-1)*4+3; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order2_3.txt
 
less ${name}_spatial_2.extract_all_crick.filter.order2.txt | awk -F "\t" '{$2=($1-1)*4+4; print $2}' > ${name}_spatial_2.extract_all_crick.filter.order2_4.txt


cat ${name}_spatial_2.extract_all_crick.filter.order2_1.txt  ${name}_spatial_2.extract_all_crick.filter.order2_2.txt  ${name}_spatial_2.extract_all_crick.filter.order2_3.txt  ${name}_spatial_2.extract_all_crick.filter.order2_4.txt | sort -k1,1n > ${name}_spatial_2.extract_all_crick.filter.order2_all.txt

rm -fr ${name}_spatial_2.extract_all_crick.filter.order2_1.txt ${name}_spatial_2.extract_all_crick.filter.order2_2.txt ${name}_spatial_2.extract_all_crick.filter.order2_3.txt ${name}_spatial_2.extract_all_crick.filter.order2_4.txt

zcat  ${Input_path1}/${name}_R1.fq.gz | less | awk 'NR==FNR{a[NR]=$0}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_crick.filter.order2_all.txt | gzip  > ${name}_spatial_Read1.extract_crick.fq.gz

rm -fr ${name}_spatial_2.extract_all_crick.filter.order2_all.txt



