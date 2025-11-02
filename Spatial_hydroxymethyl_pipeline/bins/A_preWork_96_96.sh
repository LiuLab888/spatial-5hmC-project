#!/bin/sh

usage="sh <script.sh> <data> <core> <ppn> <mem>"

if [ $# -ne 4 ]; then
        echo "$usage"
        exit
fi
Data=$1
core=$2
ppn=$3
mem=$4

source activate  /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd 

R_path=/software/biosoft/software/python/anaconda3-python3-2018/bin/Rscript
mkdir -p /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${Data}
work=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${Data}

name=${Data}_spatial
RNA_R1=/p300s/liujiang_group/tangym/Analysis_data/${Data}/${Data}_R1.fq.gz
RNA_R2=/p300s/liujiang_group/tangym/Analysis_data/${Data}/${Data}_R2.fq.gz

mkdir -p ${work}/spatial_${Data}
mkdir -p ${work}/spatial_${Data}/fastq
out_path=${work}/spatial_${Data}


##96_96 Barcode
ID_File=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation/combine_barcode.round2round1_index1_index2.v3_big_methylation.txt


#Taking about 10min

cat > spatial_${Data}.sh << EOF

#!/bin/sh
#PBS -N ${Data}_pre
#PBS -o ${work}/${Data}_pre.o
#PBS -e ${work}/${Data}_pre.e
#PBS -q ${core}
#PBS -l mem=${mem}gb,walltime=999:00:00,nodes=1:ppn=${ppn}
#HSCHED -s hschedd
#PPN limit ${ppn}

##1
sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/1_extract_barcode_methylation.sh $name $RNA_R1 $RNA_R2 $out_path/fastq

##2
zcat $out_path/fastq/${name}_2.extract.fq.gz | awk 'NR%4==2{print substr(\$1,1,22)}' > $out_path/fastq/${name}_2.extract.barcode

less $out_path/fastq/${name}_2.extract.barcode | awk -F "\t" '{\$2=NR; print \$0}' OFS="\t" > $out_path/fastq/${name}_2.extract.barcode.order

$R_path /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/2_check_barcode.R -b $ID_File -p $out_path/fastq -n ${name}_2.extract.barcode


EOF

#dsub spatial_${Data}.sh


