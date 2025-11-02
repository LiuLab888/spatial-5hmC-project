#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
from collections import defaultdict
desc = ''' '''
parser = argparse.ArgumentParser(prog = sys.argv[0] + os.linesep,description=desc, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('line_num', nargs='?', help = 'line num and fq for pick fq gzip')
parser.add_argument('fqgz', nargs='?', help = 'line num and fq for pick fq gzip')
parser.add_argument('-s', nargs='?', help = 'span for line bin', default = 1000, type = int)
parser.add_argument('-o', nargs='?', help = 'output file name', required = True)
if len(sys.argv) == 1:
    parser.print_help().__str__
    sys.exit(2)
args = parser.parse_args()



def line_parse():
    infor = defaultdict( list )
    with open( args.line_num ) as f :
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            read_num = int(line) // 4
            infor[ read_num // args.s ].append( read_num )
    return infor

def pick_fq_to_gz( reads_num ):
    with gzip.open( args.fqgz, 'rb' ) as f :
        ofh = gzip.open( args.o, 'wb')
        for i,v in enumerate( f ):
            ibin = i // args.s
            l2,l3,l4 = [ next(f) for i in range(3) ]
            if i in reads_num[ibin] and i == reads_num[ibin][0]:
                for e in [v, l2, l3, l4]:
                    ofh.write( e )
                reads_num[ibin].pop(0)
        ofh.close()

if __name__ == '__main__':
    reads_num = line_parse()
    pick_fq_to_gz( reads_num )





























