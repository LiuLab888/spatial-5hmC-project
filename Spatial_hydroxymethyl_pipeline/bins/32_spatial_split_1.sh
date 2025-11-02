#!/bin/sh
usage="sh <script.sh> <sample_name> <queue> <ppn> <mem> <thread> <part> <column>"

if [ $# -ne 7 ]; then
        echo "$usage"
        exit
fi

name=$1
queue=$2
ppn=$3
mem=$4
n=$5 # ppn:n=2:1
part=$6
m=$7


combine_barcode=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation/only_barcode_one_row.txt

Output_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}

cd ${Output_path}

work_path=${Output_path}/fastq

mkdir -p ${Output_path}/bismark/sorting/${part}

Output_path1=${Output_path}/bismark/sorting/${part}


#Taking about 3 hours for around 20000000 reads/2G of file


source activate /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd

cat > ${name}_${part}_split1.sh << EOF

#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${Output_path}/split1_${part}.err
#PBS -o ${Output_path}/split1_${part}.out
#HSCHED -s hschedd
#PPN limit ${ppn}


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

sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/spatial_split_core_1.sh ${name} \${i} ${work_path}

echo >&3

}&

done

wait

EOF

dsub ${name}_${part}_split1.sh


