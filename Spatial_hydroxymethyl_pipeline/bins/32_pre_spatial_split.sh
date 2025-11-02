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


work_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/fastq

work_path1=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/bismark/integration

cd ${work_path}

source activate  /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd/hschedd

cat > ${name}_pre_split.sh << EOF



#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${work_path}/${name}_pre_split.err
#PBS -o ${work_path}/${name}_pre_split.out
#HSCHED -s hschedd
#PPN limit ${ppn}


cd ${work_path}


gunzip -c ${name}_spatial_1.extract.fq.gz > ${name}_spatial_1.extract.txt


less ${name}_spatial_1.extract.txt | awk 'NR%4==1{print \$0}' > ${name}_spatial_1.extract_index.txt


cd ${work_path1}


samtools view -H ${name}.SE.merge.bam > header.txt


#filter

samtools view ${name}.SE.merge.bam | awk '{print \$1,\$14}' | sed 's/\.//g' > ${name}.bam.step1.txt


wc -l ${name}.bam.step1.txt > Filter.bam.summary


cat ${name}.bam.step1.txt | awk -F [":"," "] '{print \$13}' | sed s/[X,H,Z]/X/g > ${name}.bam.step2.txt


paste ${name}.bam.step1.txt ${name}.bam.step2.txt > ${name}.bam.step3.txt


cat ${name}.bam.step3.txt | awk -F ["\t"," ",":"] '{print \$1":"\$2":"\$3":"\$4":"\$5":"\$6":"\$7":"\$8":"\$9":"\$10" "\$11":"\$12":"\$14}' > ${name}.bam.step4.txt


grep -v XXX ${name}.bam.step4.txt | awk -F " " '{print \$1}' > ${name}.bam.step5.txt


wc -l ${name}.bam.step5.txt >> Filter.bam.summary


java -jar /software/biosoft/software/picard-tools-1.119/FilterSamReads.jar I=${name}.SE.merge.bam O=${name}.SE.merge.filtered.bam READ_LIST_FILE=${name}.bam.step5.txt FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT


samtools view ${name}.SE.merge.filtered.bam > ${name}.SE.merge.filtered.sam

split -l 10000000 -a 1 -d  ${name}.SE.merge.filtered.sam ${name}.SE.merge.filtered 

less ${name}.SE.merge.filtered0 | awk -F "\t" '{\$18=NR; print \$1"\t"\$18}' > ${name}_bismark_index0

less ${name}.SE.merge.filtered1 | awk -F "\t" '{\$18=NR; print \$1"\t"\$18}' > ${name}_bismark_index1



EOF

dsub ${name}_pre_split.sh



