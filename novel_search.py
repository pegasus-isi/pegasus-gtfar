#!/usr/bin/env python

import sys
from modules.NovelGene  import *
from modules.ProgressBar  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''




def novelCall(sortedLocations,prefix,key,species,VERBOSE):

  
   
    
    progressBar = ProgressBar(sys.stdout,"Searching for novel genes...",'.','Complete',1,VERBOSE)
    progressBar.increment()


    novel = NovelGene(sortedLocations,key,prefix,species)

    while novel.open: 
        novel.search()
    
    progressBar.complete()

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-s", "--chooseSpecies", default = "HUMAN", type='string', help="Change Species")
    parser.add_option("-k", "--key", default = None, type='string', help="gene key")
    parser.add_option("-p", "--prefix", default = "test", type='string', help="prefix file name")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    
    (options, args) = parser.parse_args()

    if options.chooseSpecies.upper() in ["MONKEY","RHESUS"]:          SPECIES="RHESUS"
    else:                                                           SPECIES="HUMAN"



    if len(args) == 1 and options.key!=None:
        novelCall(args[0],options.prefix,options.key,options.chooseSpecies,options.verbose)
    else:

        parser.print_help()
        print ""
        print "Example Usage: ./novel_search.py  gene_map.srt -k mykey --verbose"
        sys.exit(2)

