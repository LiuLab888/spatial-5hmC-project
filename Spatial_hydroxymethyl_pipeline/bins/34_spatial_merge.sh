#!/bin/sh
usage="sh <script> <work_path>"

if [ $# -ne 1 ]; then
	
	echo "$usage"
	
	exit

fi


work_path=$1




cat ${work_path}/part_*/methylation_matrix.txt  > ${work_path}/methylation_matrix_total.txt


cat ${work_path}/part_*/coverage_matrix.txt  > ${work_path}/coverage_matrix_total.txt


cat ${work_path}/part_*/depth_matrix.txt  > ${work_path}/depth_matrix_total.txt
