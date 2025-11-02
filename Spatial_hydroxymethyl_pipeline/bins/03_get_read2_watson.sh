#!/bin/bash
usage="sh <script.sh> <name> <Input_path> "

if [ $# -ne 2 ]; then

    echo "$usage"

    exit
fi


name=$1

Input_path=$2



cd ${Input_path}

zcat ${name}_spatial_2.extract_all_watson.fq.gz | awk 'NR%4==1{print $0}' > ${name}_spatial_2.extract_all_watson.filter.txt

less ${name}_spatial_2.extract_all_watson.filter.txt | cut -f1 -d " " | cut -f1 -d "_" > ${name}_spatial_2.extract_all_watson.filter_1_col.txt



#zcat ${name}_spatial_2.extract.fq.gz |  awk 'NR%4==1{print \$0}' >  ${name}_spatial_2.extract.filter.txt

#less ${name}_spatial_2.extract.filter.txt  | cut -f1 -d " " | cut -f1 -d "_" > ${name}_spatial_2.extract.filter_1_col.txt 



less ${name}_spatial_2.extract.filter_1_col.txt | awk 'NR==FNR{a[$1]=NR}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_watson.filter_1_col.txt | grep -v '^$' > ${name}_spatial_2.extract_all_watson.filter.order.txt


less ${name}_spatial_2.extract_all_watson.filter.order.txt | awk -F "\t" '{$2=($1-1)*4+1; print $2}' > ${name}_spatial_2.extract_all_watson.filter.order_1.txt
 
less ${name}_spatial_2.extract_all_watson.filter.order.txt | awk -F "\t" '{$2=($1-1)*4+2; print $2}' > ${name}_spatial_2.extract_all_watson.filter.order_2.txt
 
less ${name}_spatial_2.extract_all_watson.filter.order.txt | awk -F "\t" '{$2=($1-1)*4+3; print $2}' > ${name}_spatial_2.extract_all_watson.filter.order_3.txt
 
less ${name}_spatial_2.extract_all_watson.filter.order.txt | awk -F "\t" '{$2=($1-1)*4+4; print $2}' > ${name}_spatial_2.extract_all_watson.filter.order_4.txt


cat ${name}_spatial_2.extract_all_watson.filter.order_1.txt  ${name}_spatial_2.extract_all_watson.filter.order_2.txt  ${name}_spatial_2.extract_all_watson.filter.order_3.txt  ${name}_spatial_2.extract_all_watson.filter.order_4.txt | sort -k1,1n > ${name}_spatial_2.extract_all_watson.filter.order_all.txt

rm -fr ${name}_spatial_2.extract_all_watson.filter.order_1.txt ${name}_spatial_2.extract_all_watson.filter.order_2.txt ${name}_spatial_2.extract_all_watson.filter.order_3.txt ${name}_spatial_2.extract_all_watson.filter.order_4.txt

zcat  ${name}_spatial_2.extract.fq.gz | awk 'NR==FNR{a[NR]=$0}NR!=FNR{print a[$1]}' - ${name}_spatial_2.extract_all_watson.filter.order_all.txt | gzip  > ${name}_spatial_2.extract_watson.fq.gz

rm -fr ${name}_spatial_2.extract_all_watson.filter.order_all.txt


