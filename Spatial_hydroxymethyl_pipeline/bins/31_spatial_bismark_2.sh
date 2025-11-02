#!/bin/sh
usage="sh <script.sh> <sample_name> <queue> <ppn> <mem> "

if [ $# -ne 4 ]; then
 
	echo "$usage"

	exit

fi

name=$1
queue=$2
ppn=$3
mem=$4


Output_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}

cd ${Output_path} 


bismark_index=/xtdisk/liujiang_group/tangym/Ref_genomes/index/bismark_index_for_TAB_seq


bismark=/xtdisk/liujiang_group/tangym/Softwares/Bismark-0.23.0



Output_path1=${Output_path}/fastq

mkdir -p  ${Output_path}/bismark

Output_path2=${Output_path}/bismark


#Taking about 30min

#part one

source activate /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd

cat > ${name}_bismark13.sh << EOF


#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${Output_path2}/bismark13_${name}.err
#PBS -o ${Output_path2}/bismark13_${name}.out
#HSCHED -s hschedd
#PPN limit ${ppn}

## Bismark


cd ${Output_path2}

mkdir -p  integration

cd integration



samtools cat -o ${name}.SE.merge.bam ${name}_spatial_crick_unmapped_*_trimmed_bismark_bt2.bam ${name}_spatial_watson_unmapped_*_trimmed_bismark_bt2.bam  ${name}_spatial_crick_trimmed_bismark_bt2.bam ${name}_spatial_watson_trimmed_bismark_bt2.bam




EOF

dsub ${name}_bismark13.sh





