#!/usr/bin/env python


import os
import sys
#from modules.fastq.FastqFile import *
import modules.ToolSet as gtTools 
from modules.FastqFilter import FastqReads
				
def process_file(INFO,readlen,trimLengths,strictLevel,tpAdaptors,fpAdaptors):
    
    fastq = FastqFile(INFO,readlen,trimLengths,strictLevel)
    #fastq.loadAdaptor('three_prime_head',"START",'AGATCGGAAGAGCACACGT')
    #fastq.loadAdaptor('five_prime_tail',"END","CACGACGCTCTTCCGATCT")
    #fastq.loadAdaptor('three_prime_tail',"END",'GTATGCCGTCTTCTGCTTG')
    
    sys.exit()
    while fastq.nextRead():
        fastq.adaptorSearch()
        #fastq.rateQuality()
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
    parser.add_option("-t", "--trimLengths", default = '50,75', type='string', help="Trimmable Lengths")
    parser.add_option("-q", "--lowqual", default = 1, type='int', help="Low Quality Score")
    parser.add_option("-a", "--avgqual", default = 15, type='int', help="Minimum Avg Quality Score")
    parser.add_option("-p", "--prefix", default = 'foo', type='string', help="Output Filename Prefix")
    parser.add_option("-s", "--strictness", default = 1, type='int', help="Output Filename Prefix")
    #parser.add_option("-m", "--multiple", default = 25, type='int', help="Trimming Multiple")
    parser.add_option("--threePrimeAdaptors",default ="AGATCGGAAGAGCACACGT,GTATGCCGTCTTCTGCTTG",type='string',help='three prime adaptor seqs')
    parser.add_option("--fivePrimeAdaptors",default ="ACACTCTTTCCCTA,CACGACGCTCTTCCGATCA",type='string',help='five prime adaptor seqs')
    parser.add_option("--filterSeqs",default='',type='string',help="corrupt sequences")

    (options, args) = parser.parse_args()

    try:                FILE=gtTools.fileGrab(args,0)
    except IndexError:  
        parser.print_help()
        sys.exit()
    reads = FastqReads(FILE,options.readlen,options.prefix) 
    reads.addTrimLengths(options.trimLengths)    
    reads.addAdaptors(options.threePrimeAdaptors,options.fivePrimeAdaptors,options.strictness)
    
    while reads.fileOpen:
        reads.adaptorFilter()
        reads.nextRead()
    reads.printSummaryStats()
    

