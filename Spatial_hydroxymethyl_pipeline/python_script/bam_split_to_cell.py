#!/usr/bin/env python3
import os
import re
import sys
import gzip
import argparse
import subprocess as sp
import multiprocessing as mp
from collections import defaultdict
from env_name.nSys import system
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument( 'bam', nargs='?', help ='bam file')
parser.add_argument( 'cell_fqgz', nargs = '?', help = 'sample.extract_R2.fq.gz')
parser.add_argument( '-ref', nargs = '?', help = 'cell_index only_barcode_one_row.txt', default = '/mnt/server1/data3/tangym/p300s/tangym/Analysis_data/spatial_Hydroxymethylation/only_barcode_one_row.txt')
parser.add_argument( '-c', nargs = '?', help = 'cpu for calculate', default = 6, type = int)
parser.add_argument( '-e', nargs = '?', help = 'barcode each part. default 500', default = 500, type = int)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()


def barcodes_parse():
    barcodes, p = set(), re.compile(r'\s+')
    with open( args.ref ) as f :
        for line in f:
            line_arr = p.split(line.strip())
            for e in line_arr :
                e_cells = p.split(e)
                for r in e_cells:
                    if r :
                        barcodes.add(r)
    print ( 'barcode num is: ', len(barcodes), file = sys.stderr)
    return barcodes
def fq_filter_by_refrence( barcodes ):
    infor = defaultdict( str )
    infor_reads = defaultdict( int )
    line_show = 5
    with gzip.open( args.cell_fqgz ) as f :
        for read in f :
            read = read.decode().strip().strip('@').split(' ')
            read = read[0].split('_')[0]
            barcode = f.readline().decode().strip()
            f.readline()
            f.readline()
            if barcode in barcodes :
                infor_reads[barcode] += 1
                infor[ read ] = barcode
    with open( os.path.basename( args.bam ) +'.log', 'w' ) as f :
        for barcode in infor_reads:
            print ( barcode, 'reads num:', infor_reads[barcode], file = f)
    return infor

def split_bam_to_cell( infor, barcodes ):
    if 1 :
        header = system.run( 'samtools view -H {}'.format( args.bam ), shell = True )
        barcode_ofh = defaultdict(str)
        j = 0
        for i, barcode in enumerate(barcodes) :
            if not i % args.e :
                j += 1
                barcode_new_dir = system.dir( os.path.join( os.path.abspath('.'), 'part{}'.format(j) ) ).check()
            barcode_ofh[barcode] = open( os.path.join(barcode_new_dir,barcode +'.sam'), 'w' )
            print ( *header, sep = '\n', file = barcode_ofh[barcode] )
    cmd = tuple('samtools view {}'.format(args.bam).split(' ') )
    print ( cmd, file = sys.stderr )
    p = sp.Popen( cmd, stdout = sp.PIPE, stderr = sp.PIPE, bufsize = 1, text = True )
    line_show = 5
    for line in p.stdout:
        line = line.strip()
        if not line :
            continue
        line_arr = line.split(' ')
        read = line_arr[0].split('_')[0]
        if read in infor :
            barcode = infor[ read ]
            infor.pop(read)
            print ( line, file = barcode_ofh[ barcode ] )
    for barcode in barcode_ofh:
        barcode_ofh[barcode].close()
    with mp.Pool( args.c ) as p:
        p.map( sort, [ i.name for i in barcode_ofh.values() ] )
    p.close()
    p.join()
def sort( sam ):
    cmd = 'samtools view -bS {} > {}'.format( sam, sam.replace('.sam','.bam'))
    os.system(cmd)
    os.system('rm {}'.format(sam))

if __name__ == '__main__':
    #split_bam_to_cell( {}, () )
    barcodes = barcodes_parse()
    infor = fq_filter_by_refrence( barcodes )
    split_bam_to_cell( infor, barcodes )


























