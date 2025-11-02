#!/usr/bin/env python3
import os
import sys
import re
import gzip
import argparse
import subprocess as sp
import multiprocessing as mp
from collections import defaultdict
from env_name.nSys import system
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'd', nargs='?', help = 'dir for find bams calculate')
parser.add_argument( '-b', nargs='?', help = 'bismark dir', default = '/home/tangym/software/Bismark-0.23.0/')
parser.add_argument( '-bi', nargs='?', help = 'bismark index', default = '/mnt/server1/data3/tangym/xtdisk/tangym/Ref_genomes/index/bismark_index_for_TAB_seq')
parser.add_argument( '-c', nargs='?', help = 'cpu for multi', default = 6, type = int)
parser.add_argument( '-cm', action='store_true', help = 'calculate the level and merge to matrix')
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def bismark( bam ):
    prefix = bam.replace('.bam','')
    cmd = 'cd {} && deduplicate_bismark -s --bam {} 2>/dev/null 1>/dev/null'.format( args.d, bam)
    os.system(cmd)
    dup_bam = '{}.deduplicated.bam'.format(prefix)
    cmd = 'cd {} && bismark_Hydroxymethylation_extractor -s --multicore 1 --bedGraph --genome_folder {} {} 2>/dev/null 1>/dev/null'.format( args.d, args.bi, dup_bam)
    os.system(cmd)
    cmd = 'rm -fr CpG*{}* CHH*{}* CHG*{}*'.format(prefix,prefix,prefix)
    #os.system(cmd)
    cmd = 'rm -fr {}*.M-bias.txt {}*_splitting_report.txt'.format(prefix,prefix)
    #os.system(cmd)
def merge( ):
        d = args.d
        cov_gzs = [ i for i in os.listdir(d) if i.endswith('.deduplicated.bismark.cov.gz') ]
        print ('Find {} .deduplicated.bismark.cov.gz fils'.format(len(cov_gzs)))
        chroms = list(map( str, range(30)))
        chroms.extend(['X','Y','MT'])
        Hydroxymethylation_matrix = open(os.path.join(d,'Hydroxymethylation_matrix.txt'),'w')
        cov_matrix = open(os.path.join(d,'coverage_matrix.txt'),'w')
        for gz in cov_gzs:
            mC, C, cov = 0, 0, 0
            barcode = os.path.basename(gz).replace('.deduplicated.bismark.cov.gz','')
            with gzip.open(gz) as f :
                for line in f :
                    line_arr = line.strip().decode().split('\t')
                    if line_arr[0].replace('chr','') not in chroms:
                        continue
                    try :
                        line_arr[4:] = map( int, line_arr[4:] )
                    except :
                        print ( line_arr, gz )
                        exit()
                    if line_arr[4] > 0 or line_arr[5] > 0 :
                        cov += 1
                        mC += line_arr[4]
                        C += line_arr[4] + line_arr[5]
            methlyation_level = mC/C if C else 'NaN'
            print ( barcode, methlyation_level, sep = '\t', file = Hydroxymethylation_matrix)
            print ( barcode, cov, sep = '\t', file = cov_matrix)
        Hydroxymethylation_matrix.close()
        cov_matrix.close()
        bam_suff = '.bam'
        depth_matrix = open(os.path.join(d,'depth_matrix.txt'),'w')
        all_bams = [ i for i in os.listdir(d) if i.endswith(bam_suff) and len(i.split('.')) == 2 ]
        for bam in all_bams:
            depth = list(system.run('samtools view {} | wc -l'.format(bam), shell = True))[0]
            print ( bam.replace( bam_suff, ''), depth, sep = '\t', file = depth_matrix)
        depth_matrix.close()

def check_already_given_bams( ):
        d = args.d
        suffix, bam_suff = '.deduplicated.bismark.cov.gz','.bam'
        already_barcode = [ i.split('.')[0] for i in os.listdir(d) if i.endswith(suffix) and os.path.getsize(i) > 0 and len(i.split('.')) == 5 ]
        all_bams = [ i for i in os.listdir(d) if i.endswith(bam_suff) and len(i.split('.')) == 2 ]
        #print ( [ i.replace( bam_suff, '') for i in all_bams if i.replace( bam_suff, '') not in already_barcode ] )
        need = [ os.path.abspath(i) for i in all_bams if i.split('.')[0] not in already_barcode ]
        print ( 'Already: {}, all bams: {}, need calculate: {}'.format(len(already_barcode), len(all_bams), len(need)), file = sys.stderr )
        print ( need[:3] )
        return need
if __name__ == '__main__' :
    kwargs = vars( args )
    #dirs = list(map( str.strip, args.d.readlines() ) )
    #dirs = [ os.path.abspath(i) for i in dirs ]
    #kwargs.update({'dirs': dirs })
    if not args.cm:
        bams = check_already_given_bams()
        with mp.Pool( args.c ) as p:
            p.map( bismark, bams )
        p.close()
        p.join()
    else :
        merge()








