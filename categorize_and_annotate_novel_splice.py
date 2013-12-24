#!/usr/bin/env python


import sys
import os
from collections import defaultdict as dd	
from modules.Utilities import *
from modules.SpliceRecord import *
from modules.GtfFile import *
from math import fabs


def process_files(args,spliceKey,gtfFile,genomePath,readlen):

    spliceRecord = SpliceRecord(1)
    spliceRecord.addKey(spliceKey)
    spliceRecord.addGappedAlignmentFiles(args)

    gtf = GtfFile(gtfFile,genomePath,readlen)
    gtf.addSpliceRecord(spliceRecord)
    while gtf.open:
        gtf.loadGenesOnChromosome()
        gtf.addFasta(genomePath+'/'+gtf.chr+'.fa')
        gtf.printSpliceData()
        gtf.startNextChromosome()







                            

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-e", "--expr",  action = 'store_true', default = False, help="expression cnts")
    parser.add_option("-s", "--splice", action = 'store_true', default = False, help="splice cnts")
    parser.add_option("-k", "--key", action = 'store', default = None, help="phenotype key")
    parser.add_option("-n", "--norm", action = 'store', default = "MEAN", help="normalize")
    parser.add_option("-r", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-g", "--genomePath", default = 'None', type='string', help="Path to genome chr fasta files")
    parser.add_option("-f", "--gtf", default = None, type='string', help="Path to gtf file")


    (options, args) = parser.parse_args()

    if len(args)<1:
        print "WE NEED MORE FILES" 
        print args
        sys.exit()
    elif not options.key:
        print "A key is required"
        sys.exit()
    else:
        NORM=options.norm.upper()
        process_files(args,options.key,options.gtf,options.genomePath,options.readlen)
