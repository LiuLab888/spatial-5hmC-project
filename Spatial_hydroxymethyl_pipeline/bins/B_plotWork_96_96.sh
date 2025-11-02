#!/bin/sh
usage="sh <script.sh> <name> <queue> <ppn> <mem>"

if [ $# -ne 4 ]; then

	echo "${usage}"
	
	exit

fi

name=$1
queue=$2
ppn=$3
mem=$4

R_path=/software/biosoft/software/python/anaconda3-python3-2018/bin/Rscript

work_path=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}

input_file=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/Results_data/${name}/spatial_${name}/bismark/matrix

refer_file=/p300s/liujiang_group/tangym/Analysis_data/spatial_methylation/combine_barcode.round2round1_index1_index2.v3_big_methylation.txt

cd ${work_path}

#Taking about 5min

source activate hschedd

cat > spatial_plot.sh << EOF

#!/bin/sh
#PBS -q ${queue}
#PBS -l nodes=1:ppn=${ppn},mem=${mem}gb,walltime=200:00:00
#PBS -e ${work_path}/plot_${name}.err
#PBS -o ${work_path}/plot_${name}.out
#HSCHED -s hschedd
#PPN limit ${ppn}


sh /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/34_spatial_merge.sh ${work_path}/bismark/matrix


${R_path} /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/4_plot_methylation_matrix.R -p  ${work_path} -n  ${input_file}/methylation_matrix_total.txt -b  ${refer_file}


${R_path} /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/5_plot_coverage_matrix.R -p ${work_path} -n ${input_file}/coverage_matrix_total.txt  -b ${refer_file}


${R_path} /p300s/liujiang_group/tangym/Analysis_data/spatial_methylation_bismark_integration_strand/6_plot_depth_matrix.R -p ${work_path} -n ${input_file}/depth_matrix_total.txt -b  ${refer_file}



EOF

dsub spatial_plot.sh


