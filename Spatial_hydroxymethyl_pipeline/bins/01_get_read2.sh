#!/bin/bash
usage="sh <script.sh> <name> <queue> <ppn> <mem> "

if [ $# -ne 4 ]; then

    echo "$usage"

    exit
fi


name=$1

queue=$2

ppn=$3

mem=$4



Input_path1=/p300s/liujiang_group/tangym/Analysis_data/${name}

Input_path2=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/fastq



source activate hschedd

cd ${Input_path2}

cat > ${name}_get_read2.sh << EOF


#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${Input_path2}/get_read2_${name}.err
#PBS -o ${Input_path2}/get_read2_${name}.out
#HSCHED -s hschedd
#PPN limit ${ppn}


cd ${Input_path2}

zcat ${name}_spatial_2.extract.fq.gz| less | awk 'NR%4==1{print \$0}' > ${name}_spatial_2.extract.filter.txt

less ${name}_spatial_2.extract.filter.txt | cut -f1 -d " " | cut -f1 -d "_" > ${name}_spatial_2.extract.filter_1_col.txt


zcat ${Input_path1}/${name}_R2.fq.gz | less | awk 'NR%4==1{print \$0}' >  ${name}_R2.filter.txt

less ${name}_R2.filter.txt  | cut -f1 -d " " > ${name}_R2.filter_1_col.txt 



less ${name}_R2.filter_1_col.txt | awk 'NR==FNR{a[\$1]=NR}NR!=FNR{print a[\$1]}' - ${name}_spatial_2.extract.filter_1_col.txt | grep -v '^$' > ${name}_spatial_2.extract.filter.order.txt


less ${name}_spatial_2.extract.filter.order.txt | awk -F "\t" '{\$2=(\$1-1)*4+1; print \$2}' > ${name}_spatial_2.extract.filter.order_1.txt
 
less ${name}_spatial_2.extract.filter.order.txt | awk -F "\t" '{\$2=(\$1-1)*4+2; print \$2}' > ${name}_spatial_2.extract.filter.order_2.txt
 
less ${name}_spatial_2.extract.filter.order.txt | awk -F "\t" '{\$2=(\$1-1)*4+3; print \$2}' > ${name}_spatial_2.extract.filter.order_3.txt
 
less ${name}_spatial_2.extract.filter.order.txt | awk -F "\t" '{\$2=(\$1-1)*4+4; print \$2}' > ${name}_spatial_2.extract.filter.order_4.txt


cat ${name}_spatial_2.extract.filter.order_1.txt  ${name}_spatial_2.extract.filter.order_2.txt  ${name}_spatial_2.extract.filter.order_3.txt  ${name}_spatial_2.extract.filter.order_4.txt | sort -k1,1n > ${name}_spatial_2.extract.filter.order_all.txt

rm -fr ${name}_spatial_2.extract.filter.order_1.txt ${name}_spatial_2.extract.filter.order_2.txt ${name}_spatial_2.extract.filter.order_3.txt ${name}_spatial_2.extract.filter.order_4.txt

zcat  ${Input_path1}/${name}_R2.fq.gz | less | awk 'NR==FNR{a[NR]=\$0}NR!=FNR{print a[\$1]}' - ${name}_spatial_2.extract.filter.order_all.txt | gzip  > ${name}_spatial_2.extract_all.fq.gz

rm -fr ${name}_spatial_2.extract.filter.order_all.txt

EOF

dsub ${name}_get_read2.sh


