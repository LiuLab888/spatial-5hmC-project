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


Out_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/fastq

software_path=/asnas/liujiang_group/yuhao/Software/anaconda3/bin

source activate hschedd

cat > ${name}_get_strand_crick.sh << EOF

#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${Out_path}/get_strand_crick_${name}.err
#PBS -o ${Out_path}/get_strand_crick_${name}.out
#HSCHED -s hschedd
#PPN limit ${ppn}

cd ${Out_path}

${software_path}/umi_tools extract --extract-method=regex \
--bc-pattern="^(?P<discard_1>GGTGTAGTGGGTTTGGAGG){s<=3}(?P<cell_1>.{11})(?P<discard_2>..TT.T...T.....T.T.T..T...T...){s<=4}(?P<cell_2>.{11})(?P<discard_3>TTT.....T..TT....T...T...TT...){s<=4}(?P<discard_4>AGATGTGT){s<=2}(?P<umi_1>.)(?P<discard_5>.*)$" \
-I ${name}_spatial_2.extract_all.fq.gz -S ${name}_spatial_2.extract_all_crick.fq.gz \
-L ${name}_extract_crick.log

sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/separate_strand_new/03_get_read2_crick.sh  ${name}  /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/fastq

sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/separate_strand_new/04_get_crick.sh  ${name}  /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/fastq  /p300s/liujiang_group/tangym/Analysis_data/${name}


EOF

dsub ${name}_get_strand_crick.sh


