#!/usr/bin/env python


import sys
import os
from collections import defaultdict as dd	
from modules.Utilities import *
from modules.SpliceRecord import *
from modules.GtfFile import *
from math import fabs








def process_files(myfiles):
    splice_cnt = 0; uniq_cnt =0; total =0
    for f in myfiles:
        INIT=True
        try:
            handle = open(f)
        except IOError:
            errorQuit("Invalid File "+f)

        for line in handle:
            if line[0]=="@": continue
            line=line.split()
            if INIT:
                UNIQ=True
                ID=line[0];
                if len(line[5].split("M"))>2: SPLICED=True
                else:                       SPLICED=False
                INIT=False
                continue 
            else:
                if line[0] == ID:
                    UNIQ=False
                    if len(line[5].split("M"))>2: SPLICED=True
                    else:                       SPLICED=False
                else:
                    total+=1
                    if UNIQ:
                        uniq_cnt+=1
                        if SPLICED: splice_cnt+=1
                    UNIQ=True; ID=line[0]
                    if len(line[5].split("M"))>2: SPLICED=True
                    else:                       SPLICED=False
    
    print "QC-Statistics: GENOMIC MAPPING"
    print "Total-Reads[total,uniq/repetitive/ambiguous]", total, uniq_cnt, total-uniq_cnt, 0
    print "Uniq-Reads[spliced/unspliced]", uniq_cnt, splice_cnt 



if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    #parser.add_option("-k", "--key", action = 'store', default = None, help="phenotype key")
    #parser.add_option("-u", "--ungapped", default = None, type='string', help="Path to ungapped genomic file")
    #parser.add_option("-g", "--gapped", default = None, type='string', help="Path to gapped genomic file")


    (options, args) = parser.parse_args()

    process_files(args)
