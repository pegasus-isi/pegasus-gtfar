#!/usr/bin/env python


import os
import sys
#from modules.fastq.FastqFile import *
from modules.FastqFile import *
				
def process_file(INFO,readlen,lowQual,avgQual,minTrim,outname,multiple):
    
    fastq = FastqFile(INFO,readlen,lowQual,avgQual,minTrim,multiple)
    fastq.loadAdaptor('three_prime_head',"START",'AGATCGGAAGAGCACACGT')
    fastq.loadAdaptor('five_prime_tail',"END","CACGACGCTCTTCCGATCT")
    fastq.loadAdaptor('three_prime_tail',"END",'GTATGCCGTCTTCTGCTTG')
    
    while fastq.nextRead():
        fastq.adaptorSearch()
        fastq.rateQuality()
        if fastq.filterPass():
            fastq.writeFastq(outname)
        else:
            fastq.tryTrimFastq(outname,minTrim)
        fastq.recordStats()
    fastq.printStats()
                
                
            


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Read Length")
    parser.add_option("-t", "--trimlen", default = 50, type='int', help="Minimum Trim Length")
    parser.add_option("-q", "--lowqual", default = 1, type='int', help="Low Quality Score")
    parser.add_option("-a", "--avgqual", default = 15, type='int', help="Minimum Avg Quality Score")
    parser.add_option("-p", "--prefix", default = 'foo', type='string', help="Output Filename Prefix")
    parser.add_option("-m", "--multiple", default = 25, type='int', help="Trimming Multiple")

    (options, args) = parser.parse_args()

    if len(args) == 1:
        FILE=open(args[0])
    else:
        FILE=sys.stdin

    process_file(FILE,options.readlen,options.lowqual,options.avgqual,options.trimlen,options.prefix,options.multiple)

