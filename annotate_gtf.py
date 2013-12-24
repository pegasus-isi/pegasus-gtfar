#!/usr/bin/env python

import sys
import os
from modules.GtfFile  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''


def runthrough(fName,genomePath,filters,prefix,readlen,mutationList):
   
    gtf = GtfFile(fName,prefix,readlen,filters,mutationList,printKEY=True)
    while gtf.open:
        gtf.loadGenesOnChromosome()
        gtf.addFasta(genomePath+'/'+gtf.chr+'.fa') 
        gtf.printGenesOnChromosome()
        gtf.startNextChromosome()
    

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-l", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mygtf', type='string', help="Output Filename Prefix")
    parser.add_option("-c", "--chromosomes", default = None, type='string', help="Path to genome autosome fasta files")
    parser.add_option("-f", "--filterseqs", default = None, type='string', help="Path to filter seqs (ribosomal/mitochondrial)")
    parser.add_option("-m", "--mutationList", default = None, type='string', help="Path mutation list")

    (options, args) = parser.parse_args()

if len(args)==1 and options.chromosomes!=None and options.filterseqs!=None:
    if options.chromosomes[-1]=="/": options.chromosomes=options.chromosomes[0:len(options.chromosomes)-1]
    runthrough(args[0],options.chromosomes,options.filterseqs,options.prefix,options.readlen,options.mutationList)
else:
    print "ANNOTATE AND CREATE SEQS/KEY"
    print "USAGE: ./annotate_gtf.py file.gtf -c chromosome_path -f filter_seqs -m mitochondria -p PREFIXNAME -r READLENGTH -m Mutationlist" 

