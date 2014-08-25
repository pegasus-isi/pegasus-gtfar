#!/usr/bin/env python

import sys
import os
from modules.GtfFile  import *
from collections import defaultdict as dd
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''


def runthrough(fName,genomePath,parameters,gene_key,filters):
    gtf = GtfFile(fName,None,None,filters)
    gtf.add_simulation_parameters(parameters)
    gtf.add_gene_key(gene_key)
    while gtf.open:
        gtf.loadGenesOnChromosome()
        gtf.addFasta(genomePath+'/'+gtf.chr+'.fa') 
        gtf.simulateGeneReadsOnChromosome()
        gtf.startNextChromosome()
    

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-l", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mygtf', type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genome", default = None, type='string', help="Path to genome autosome fasta files")
    parser.add_option("-c", "--config", default = None, type='string', help="Path to configuration file")
    parser.add_option("-k", "--key", default = None, type='string', help="Path to key file")
    parser.add_option("-f", "--filters", default = None, type='string', help="Path to filter file")

    (options, args) = parser.parse_args()




if len(args)==1 and options.genome!=None:
    if options.genome[-1]=="/": options.genome=options.genome[0:len(options.genome)-1]
    runthrough(args[0],options.genome,options.config,options.key,options.filters)
else:
    print "ANNOTATE AND CREATE SEQS/KEY"
    print "USAGE: ./annotate_gtf.py file.gtf -g chromosome_path -c file.config -k options.key" 

