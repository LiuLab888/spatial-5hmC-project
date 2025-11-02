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


source activate /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd

out_path=/p300s/liujiang_group/tangym/Analysis_data/${name}

RNA_R1=${out_path}/${name}_R1.fq.gz

RNA_R2=${out_path}/${name}_R2.fq.gz

software_path=/asnas/liujiang_group/yuhao/Software/anaconda3/bin

cd ${out_path}

cat > ${name}_whitelist.sh << EOF


#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${out_path}/whitelist_${name}.err
#PBS -o ${out_path}/whitelist_${name}.out
#HSCHED -s hschedd
#PPN limit ${ppn}


${software_path}/umi_tools whitelist --extract-method=regex \
--bc-pattern="^(?P<discard_1>.{19})(?P<cell_1>.{11})(?P<discard_2>.{30})(?P<cell_2>.{11})(?P<umi_1>.)" \
-I ${RNA_R2} \
-S ${out_path}/${name}_extract_R2.whitelist \
-L ${out_path}/${name}_extract_whitelist.log \
--error-correct-threshold=2 \
--ed-above-threshold=discard \
--set-cell-number=9216


EOF

dsub ${name}_whitelist.sh

