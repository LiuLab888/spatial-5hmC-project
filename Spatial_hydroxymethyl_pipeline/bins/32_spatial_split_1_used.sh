#!/bin/sh
usage="sh <script.sh> <rootdir> <sample_name> <part> <column>"

if [ $# -ne 4 ]; then
        echo "$usage"
        exit
fi
root=$1
name=$2
part=$3
m=$4
n=3
#sh ${script} ${name} ${queue} 20 20 ${n} part_1 1-1000
combine_barcode=/mnt/server1/data3/tangym/p300s/tangym/Analysis_data/spatial_methylation/only_barcode_one_row.txt
Output_path=$root/Results_data/${name}/spatial_${name}
bismark=/home/tangym/software/Bismark-0.23.0
Output_path=$root
sorting_path=${Output_path}/bismark/sorting/${part}
work_path=${Output_path}/fastq
fastq_dir=$work_path
mkdir -p $sorting_path
integration_dir=${Output_path}/bismark/integration


cd $sorting_path
Output_path1=${Output_path}/bismark/sorting/${part}


#Taking about 3 hours for around 20000000 reads/2G of file


cat > ${name}_${part}_split1_2.sh << EOF
#!/usr/bin/env bash
#make a tube

[ -e /tmp/fd1 ] || mkfifo /tmp/fd1

#make a signal

exec 3<>/tmp/fd1

#delete the tube

rm -f -r  /tmp/fd1



for ((k=1;k<=${n};k++))

do

echo >&3

done


cd ${Output_path1}

for i in `cut -f${m}  ${combine_barcode}`

do

read -u3

{

index=\$i
less  $fastq_dir/${name}.extract_R2.fq.barcode.order | awk -F "\\t" -v barcode=\${index} '{if(\$1==barcode) print \$2}' > \${index}.txt

less  $fastq_dir/${name}_spatial_1.extract_index.txt  | awk 'NR==FNR{a[FNR]=\$0}NR!=FNR{print a[\$1]}' - \${index}.txt  >  \${index}_fastq.txt

less \${index}_fastq.txt | awk -F " " '{\$3=substr(\$1,2,length(\$1)-2);print \$0}' OFS=" " | awk -F " " '{\$4=\$3\$2; print \$4}' > \${index}_fastq_subtract_1.txt
rm -fr \${index}_fastq.txt \${index}.txt
touch \${index}.sam
less $integration_dir/${name}_bismark_index0 | awk 'NR==FNR{a[\$1]=\$2}NR!=FNR{print a[\$1]}' -  \${index}_fastq_subtract_1.txt | grep -v '^$' > \${index}_bismark_index.order0
less $integration_dir/${name}.SE.merge.filtered0 | awk 'NR==FNR{a[NR]=\$0}NR!=FNR{print a[\$1]}' - \${index}_bismark_index.order0 >> \${index}.sam
less  \${index}.sam | grep -v '^$' >  \${index}_final.sam

cat ${integration_dir}/header.txt \${index}_final.sam > \${index}_header.sam 

samtools view -bS \${index}_header.sam > \${index}.bam

$bismark/deduplicate_bismark -s --bam \${i}.bam 
$bismark/bismark_methylation_extractor -s   --multicore 1  --bedGraph   --genome_folder  $bismark_index  \${i}.deduplicated.bam

j=\$i
rm -fr CpG*\${i}*  CHH*\${i}*  CHG*\${i}*

rm -fr \${i}*.M-bias.txt  \${i}*_splitting_report.txt


samtools view \${i}.deduplicated.bam | less | awk -v barcode=\${i} 'END{print barcode"\\t"NR}' >> depth_matrix.txt


less \${j}.deduplicated.bismark.cov.gz | awk -F "\\t"  '{if(\$1=="1" || \$1=="2" || \$1=="3" || \$1=="4" || \$1=="5" || \$1=="6" || \$1=="7" || \$1=="8" || \$1=="9" || \$1=="10" || \$1=="11" || \$1=="12" || \$1=="13" || \$1=="14" || \$1=="15" || \$1=="16" || \$1=="17" || \$1=="18" || \$1=="19" || \$1=="X" || \$1=="Y" || \$1=="MT") print \$0}' OFS="\\t" | awk -F "\\t" -v barcode=\${j} 'BEGIN{C=0;mC=0}{if(\$5>0 || \$6>0){mC=mC+\$5;C=C+\$5+\$6}}END{print barcode"\\t"mC/C}' >> methylation_matrix.txt


less \${j}.deduplicated.bismark.cov.gz | awk -F "\\t"  '{if(\$1=="1" || \$1=="2" || \$1=="3" || \$1=="4" || \$1=="5" || \$1=="6" || \$1=="7" || \$1=="8" || \$1=="9" || \$1=="10" || \$1=="11" || \$1=="12" || \$1=="13" || \$1=="14" || \$1=="15" || \$1=="16" || \$1=="17" || \$1=="18" || \$1=="19" || \$1=="X" || \$1=="Y" || \$1=="MT") print \$0}' OFS="\\t" | awk -F "\\t" -v barcode=\${j} 'BEGIN{C=0}{if(\$5>0 || \$6>0){C=C+1}}END{print barcode"\\t"C}' >> coverage_matrix.txt




echo >&3

}&

done

wait

EOF


