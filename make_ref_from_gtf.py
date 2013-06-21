#!/usr/bin/env python

import sys
from gtfar_modules.gtf.GtfGene  import *
from gtfar_modules.gtf.GtfFile  import *

'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



def runthrough(fName,genomePath,prefix,readlen=100):
   
    
    gtf = GtfFile(fName,prefix)
    while gtf.open:
        gtf.loadChromosomeGenes()
        gtf.addFasta(genomePath+'/'+gtf.chr+'.fa')
        gtf.printSeqs()
        gtf.startNewChromosome()
    gtf.pickleOut

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'foo', type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genomePath", default = None, type='string', help="Path to genome chr fasta files")

    (options, args) = parser.parse_args()

    runthrough(args[0],options.genomePath,options.prefix,options.readlen)


