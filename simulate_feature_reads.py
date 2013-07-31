#!/usr/bin/env python

import sys
import os
#from modules.gtf.GtfGene  import *
#from modules.gtf.GtfFile  import *

from modules.GtfFile import *
from modules.GtfGene import *

'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



def runthrough(fName,genomePath,prefix,readlen=100):
   
    
    gtf = GtfFile(fName,prefix,readlen,printKEY=False)
    while gtf.open:
        gtf.loadGenesOnChromosome()
        gtf.addFasta(genomePath+'/'+gtf.chr+'.fa') 
        gtf.uniquifySeqs(SILENT=True)
        gtf.simulateReads() 
        gtf.startNextChromosome()
    

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mygtf', type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genomePath", default = None, type='string', help="Path to genome chr fasta files")

    (options, args) = parser.parse_args()

if len(args)==1 and options.genomePath!=None:

    runthrough(args[0],options.genomePath,options.prefix,options.readlen)
else:
    print "ANNOTATE AND CREATE SEQS/KEY"
    print "USAGE: ./simulate_feature_reads.py file.gtf -g GENOMEPATH -p PREFIXNAME -r READLENGTH" 

