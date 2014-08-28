#!/usr/bin/env python

import sys
import os
from modules.GtfFile  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



   






    

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-l", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mygtf', type='string', help="Output Filename Prefix")
    parser.add_option("-c", "--chromosomes", default = None, type='string', help="Path to genome autosome fasta files")
    parser.add_option("-f", "--filterType", default = "HUMAN", type='string', help="Organism for filters")
    parser.add_option("-g", "--gapcands", default = None, type='string', help="File which lists possible gaps for annotation")
    parser.add_option("-m", "--mutationList", default = None, type='string', help="Path mutation list")
    parser.add_option("-s", "--snpcands", default = None, type='string', help="Path snp cands")

    (options, args) = parser.parse_args()

if len(args) == 0:
    sys.stderr.write("A gtf is required\n")
    sys.stderr.write("EXAMPLE: ./create_RDD_candidate_blocks.py file.gtf -c chromosome_path  -p PREFIXNAME -l READLENGTH -m Mutationlist\n") 
    parser.print_help()
elif options.chromosomes == None:
    sys.stderr.write("A path to chromosome files is required\n")
    sys.stderr.write("EXAMPLE: ./create_RDD_candidate_blocks.py file.gtf -c chromosome_path  -p PREFIXNAME -l READLENGTH -m Mutationlist\n") 
    parser.print_help()
elif options.snpcands == None and options.gapcands == None:
    sys.stderr.write("A path to candidate locations is is required\n")
    sys.stderr.write("EXAMPLE: ./create_RDD_candidate_blocks.py file.gtf -c chromosome_path  -p PREFIXNAME -l READLENGTH -m Mutationlist\n") 
    parser.print_help()
else:
    if options.chromosomes[-1]=="/": options.chromosomes=options.chromosomes[0:len(options.chromosomes)-1]
    if options.snpcands != None:
        gtf = GtfFile(args[0],options.prefix,options.readlen,options.filterType)
        gtf.addCandidates(options.snpcands,"SNPS")
        while gtf.open:
            gtf.loadGenesOnChromosome()
            gtf.addFasta(options.chromosomes+'/'+gtf.chr+'.fa') 
            if gtf.SNPCANDS:
                gtf.printSnpCandsOnChromosome()
            #gtf.printGenesOnChromosome()
            gtf.startNextChromosome()



