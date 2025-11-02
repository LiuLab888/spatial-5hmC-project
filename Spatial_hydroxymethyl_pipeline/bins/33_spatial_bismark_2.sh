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
n=$5 #ppn:n=2:1
part=$6
m=$7

combine_barcode=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation/only_barcode_one_row.txt


Output_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}

cd ${Output_path} 

bismark_index=/xtdisk/liujiang_group/tangym/Ref_genomes/index/bismark_index_for_TAB_seq


bismark=/xtdisk/liujiang_group/tangym/Softwares/Bismark-0.23.0


Output_path1=${Output_path}/bismark/sorting


mkdir -p ${Output_path}/bismark/matrix


Output_path2=${Output_path}/bismark/matrix



#Taking about 1.5  hour

source activate /xtdisk/liujiang_group/tangym/anaconda3/envs/hschedd

cat > ${name}_${part}_bismark2.sh << EOF

#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${Output_path}/bismark2_${part}_${name}.err
#PBS -o ${Output_path}/bismark2_${part}_${name}.out
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



## Bismark

cd ${Output_path2}

mkdir -p  ${part}

cd ${part}

touch methylation_matrix.txt

touch coverage_matrix.txt

touch depth_matrix.txt

for j  in `cut -f${m}  ${combine_barcode}`

do

read -u3

{


sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/spatial_bismark_core_2.sh \${j} ${Output_path1}/${part} ${Output_path2}/${part} $bismark $bismark_index


echo >&3

}&

done

wait


EOF

dsub ${name}_${part}_bismark2.sh


