#!/bin/sh

usage="sh <script.sh> <sample_name> <queue> "

if [ $# -ne 2 ]; then

        echo "$usage"

        exit

fi


name=$1

queue=$2

n=10


script=32_spatial_split_3.sh


sh ${script} ${name} ${queue}  20  20 ${n} part_1 1-1000  

sh ${script} ${name} ${queue}  20  20 ${n} part_2 1001-2000  

sh ${script} ${name} ${queue}  20  20 ${n} part_3 2001-3000  

sh ${script} ${name} ${queue}  20  20 ${n} part_4 3001-4000  

sh ${script} ${name} ${queue}  20  20 ${n} part_5 4001-5000  

sh ${script} ${name} ${queue}  20  20 ${n} part_6 5001-6000  

sh ${script} ${name} ${queue}  20  20 ${n} part_7 6001-7000  

sh ${script} ${name} ${queue}  20  20 ${n} part_8 7001-8000  

sh ${script} ${name} ${queue}  20  20 ${n} part_9 8001-9216 
  

