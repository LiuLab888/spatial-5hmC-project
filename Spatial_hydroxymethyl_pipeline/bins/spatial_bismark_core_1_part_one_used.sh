#!/bin/bash
help()
{
    echo "Usage: $0 [ -g | --fqgz ]
               [ -i | --index ]
               [ -p | --prefix ]
               [ -b | --bismark ]
               [ -c | --cpu ]
               [ -f | --odir1 ]
               [ -s | --odir1 ]
               [ -h | --help  ]"
    exit 2
}
SHORT=g:,i:,p:,b:,c:,f:,s:,h:,
LONG=fqgz:,index:,prefix:,bismark:,cpu:,odir1:,odir2:,help
OPTS=$(getopt -a -n weather --options $SHORT --longoptions $LONG -- "$@")
if [ "$#" -eq 0 ]; then
  help
fi
eval set -- "$OPTS"
while :
    do
        case "$1" in
            -g | -fqgz )
            fqgz="$2"
            shift 2
            ;;
            -i | --index )
            index="$2"
            shift 2
            ;;
            -p | --prefix )
            prefix="$2"
            shift 2
            ;;
            -b | --bismark )
            bismark="$2"
            shift 2
            ;;
            -c | --cpu )
            cpu="$2"
            shift 2
            ;;
            -f | --odir1 )
            od1="$2"
            shift 2
            ;;
            -s | --odir2 )
            od2="$2"
            shift 2
            ;;
            -h | --help)
            help
            ;;
            --)
            shift;
            break
            ;;
            *)
            echo "Unexpected option: $1"
            help
            ;;
        esac
    done

j=`basename $prefix`
prefix=`basename $prefix`
#echo cd $od2
if [[ $fqgz == *"crick"* ]]; then
    bismark="${bismark}/bismark --pbat"
    echo trim_galore --cores $cpu -a CTATCTCTTATACACATCT   --clip_R1 43  --basename `basename ${prefix}` $fqgz
else 
    bismark="${bismark}/bismark"
    echo trim_galore --cores $cpu -a CTGTCTCTTATACACATCT   --clip_R1 34  --three_prime_clip_R1 9   --basename `basename ${prefix}` $fqgz
fi
echo ${bismark} --unmapped --dovetail --pbat -p $cpu $index ${prefix}_trimmed.fq.gz
#echo cd $od2
#part one
#1
echo trim_galore --cores $cpu  --nextera --three_prime_clip_R1 5   --basename ${prefix}_unmapped_1  ${prefix}_trimmed.fq.gz_unmapped_reads.fq.gz
echo mv ${prefix}_unmapped_1_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${prefix}_unmapped_1_trimmed.fq.gz
echo $bismark  --unmapped --dovetail -p $cpu $index  ${prefix}_unmapped_1_trimmed.fq.gz 
#2
echo trim_galore --cores $cpu  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_2  ${j}_unmapped_1_trimmed.fq.gz_unmapped_reads.fq.gz

echo mv ${j}_unmapped_2_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_2_trimmed.fq.gz
echo $bismark  --unmapped --dovetail -p $cpu $index  ${j}_unmapped_2_trimmed.fq.gz

#3

echo trim_galore --cores $cpu  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_3  ${j}_unmapped_2_trimmed.fq.gz_unmapped_reads.fq.gz

echo mv ${j}_unmapped_3_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_3_trimmed.fq.gz

echo $bismark  --unmapped --dovetail -p $cpu $index  ${j}_unmapped_3_trimmed.fq.gz


#4

echo trim_galore --cores $cpu  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_4 ${j}_unmapped_3_trimmed.fq.gz_unmapped_reads.fq.gz

echo mv ${j}_unmapped_4_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_4_trimmed.fq.gz


echo $bismark  --unmapped --dovetail -p $cpu $index  ${j}_unmapped_4_trimmed.fq.gz


#5

echo trim_galore --cores $cpu  --nextera --three_prime_clip_R1 5   --basename ${j}_unmapped_5  ${j}_unmapped_4_trimmed.fq.gz_unmapped_reads.fq.gz

echo mv ${j}_unmapped_5_trimmed.fq.gz_unmapped_reads_trimmed.fq.gz ${j}_unmapped_5_trimmed.fq.gz


echo $bismark  --unmapped --dovetail -p $cpu $index  ${j}_unmapped_5_trimmed.fq.gz



