#!/bin/sh
usage="sh <script.sh> <sample_name> integration_dir"

if [ $# -ne 2 ]; then

        echo "$usage"

        exit

fi

name=$1
integration_dir=$2
work_path=`readlink -f ../../fastq`
work_path1=$integration_dir

cat > ${name}_pre_split.sh << EOF
cd ${work_path}
gunzip -c ${name}.extract_R1.fq.gz > ${name}_spatial_1.extract.txt


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


java -jar /home/tangym/software/picard.jar FilterSamReads I=${name}.SE.merge.bam O=${name}.SE.merge.filtered.bam READ_LIST_FILE=${name}.bam.step5.txt FILTER=includeReadList VALIDATION_STRINGENCY=LENIENT


#samtools view ${name}.SE.merge.filtered.bam > ${name}.SE.merge.filtered.sam

#split -l 10000000 -a 1 -d  ${name}.SE.merge.filtered.sam ${name}.SE.merge.filtered 

#split_sam_and_index.py ${name}.SE.merge.filtered.sam



#fls='${name}.SE.merge.filtered*'
#for i in "\${!fls[@]}"; do 
#    less "\${fls[\$i]}" | awk -F "\t" '{\$18=NR; print \$1"\t"\$18}' > ${name}_bismark_index\$i
#    less ${name}.SE.merge.filtered1 | awk -F "\t" '{\$18=NR; print \$1"\t"\$18}' > ${name}_bismark_index1

EOF




