#!/usr/bin/env python

import sys
import os
print "HI"
os.getcwd()
from modules.GtfFile  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''


def runthrough(fName,genomePath,filters,prefix,readlen,gapcands,mutationList):
   
    gtf = GtfFile(fName,prefix,readlen,filters)
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
    parser.add_option("-f", "--filterseqs", default = "HUMAN", type='string', help="Organism for filters")
    parser.add_option("-g", "--gapcands", default = None, type='string', help="File which lists possible gaps for annotation")
    parser.add_option("-m", "--mutationList", default = None, type='string', help="Path mutation list")
    parser.add_option("-s", "--snpcands", default = None, type='string', help="Path snp cands")

    (options, args) = parser.parse_args()

if len(args) == 0:
    sys.stderr.write("A gtf is required\n")
    sys.stderr.write("EXAMPLE: ./annotate_gtf.py file.gtf -c chromosome_path -f filter_seqs -p PREFIXNAME -l READLENGTH -m Mutationlist\n") 
    parser.print_help()
elif options.chromosomes == None:
    sys.stderr.write("A path to chromosome files is required\n")
    sys.stderr.write("EXAMPLE: ./annotate_gtf.py file.gtf -c chromosome_path -f filter_seqs -p PREFIXNAME -l READLENGTH -m Mutationlist\n") 
    parser.print_help()
else:
    if options.chromosomes[-1]=="/": options.chromosomes=options.chromosomes[0:len(options.chromosomes)-1]
    if options.snpcands != None:
        print "COOL"
    elif options.snpcands == None:
        runthrough(args[0],options.chromosomes,options.filterseqs,options.prefix,options.readlen,options.gapcands,options.mutationList)
    else:
        print "ANNOTATE AND CREATE SEQS/KEY"
        print "USAGE: ./annotate_gtf.py file.gtf -c chromosome_path -f filter_seqs -p PREFIXNAME -l READLENGTH -m Mutationlist" 

