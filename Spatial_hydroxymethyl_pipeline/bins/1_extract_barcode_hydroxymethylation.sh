#!/bin/sh
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

name=$1 

RNA_R1=$2 

RNA_R2=$3 

out_path=$4 

#mkdir -p ${out_path}
## 1 ##
#umi_tools extract --extract-method=regex \
#--bc-pattern2="^(?P<umi_1>.)(?P<discard_1>GGTGTAGTGGGTTTGGAGG){s<=1}.{11}(?P<discard_2>AT.{2}A.GTG.TTGAG.G.G.TG.ATA.TTG){e<=1}.{11}(?P<discard_3>.{3}ATGAT.GT.{2}GATG.AGT.GTG.{2}ATG){s<=2,i<=1,d<=1}.{10}(?P<discard_4>.*)$" \
#-I ${RNA_R1} -S ${out_path}/${name}_1.extract.fq.gz --read2-in=${RNA_R2} --read2-out=${out_path}/${name}_2.extract.fq.gz \
#-L ${out_path}/extract.log

## 2 ##
#umi_tools whitelist --extract-method=regex \
#--bc-pattern="^(?P<discard_1>.{11})(?P<umi_1>.{8})(?P<cell_1>.{11})(?P<discard_2>.{30})(?P<cell_2>.{11})" \
#-I /p300s/liujiang_group/tangym/Analysis_data/mE105sp5H_221208/mE105sp5H_221208_R2.fq.gz \
#-S /p300s/liujiang_group/tangym/Analysis_data/mE105sp5H_221208/mE105sp5H_221208_extract_R2.whitelist \
#-L /p300s/liujiang_group/tangym/Analysis_data/mE105sp5H_221208/extract_whitelist.log \
#--error-correct-threshold=2 \
#--ed-above-threshold=discard \
#--set-cell-number=9216

## 3 ##
show umi_tools extract --extract-method=regex \
--bc-pattern2="^(?P<discard_1>GGTGTAGTGGGTTTGGAGG){s<=3}.{11}(?P<discard_2>.{30}).{11}(?P<discard_3>.{30})(?P<umi_1>.)(?P<discard_4>.*)$" \
-I ${RNA_R1} -S ${out_path}/${name}_1.extract.fq.gz --read2-in=${RNA_R2} --read2-out=${out_path}/${name}_2.extract.fq.gz \
-L ${out_path}/extract.log




