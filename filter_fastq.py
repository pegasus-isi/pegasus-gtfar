#!/usr/bin/env python


import os
import sys
from gtfar_modules.fastq.FastqFile import *
				
def process_file(INFO,readlen,lowQual,avgQual,minTrim,outname):
    
    fastq = FastqFile(INFO,readlen,lowQual,avgQual,minTrim)
    fastq.loadAdaptor('three prime head',"START",'AGATCGGAAGAGCACACGT')
    fastq.loadAdaptor('five prime tail',"END","CACGACGCTCTTCCGATCT")
    fastq.loadAdaptor('three prime tail',"END",'GTATGCCGTCTTCTGCTTG')
    
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

    (options, args) = parser.parse_args()

    process_file(args[0],options.readlen,options.lowqual,options.avgqual,options.trimlen,options.prefix)

