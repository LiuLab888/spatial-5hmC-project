#!/bin/sh
usage="sh <script.sh> <root_dir> <sample_name> <part> <column> <split>"

if [ $# -ne 5 ]; then
        echo "$usage"
        exit
fi

root=$1
name=$2
part=$3
m=$4
k=$5
n=3
#sh ${script} ${name} ${queue}  20  20 ${n} part_1 1-1000 ${k}

combine_barcode=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation/only_barcode_one_row.txt

Output_path=$root

cd ${Output_path}

work_path=${Output_path}/bismark/integration

Output_path1=${Output_path}/bismark/sorting


#Taking about 5 hours


cat > ${name}_${part}_${k}_split2.sh << EOF

#!/bin/sh

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



mkdir -p  ${Output_path1}/${part}

cd ${Output_path1}/${part}


for i in `cut -f${m}  ${combine_barcode}`

do

read -u3

{

index=\$i
touch \${index}.sam
less ${work_path}/${name}_bismark_index${k} | awk 'NR==FNR{a[\$1]=\$2}NR!=FNR{print a[\$1]}' -  \${index}_fastq_subtract_1.txt | grep -v '^$' > \${index}_bismark_index.order${k}
less ${work_path}/${name}.SE.merge.filtered${k} | awk 'NR==FNR{a[NR]=\$0}NR!=FNR{print a[\$1]}' - \${index}_bismark_index.order${k} >> \${index}.sam

#sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/spatial_split_core_2.sh ${name} \${i} ${work_path} ${k}

echo >&3

}&

done

wait

EOF




