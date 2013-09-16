#!/usr/bin/env python

import sys
from modules.MutationRecord  import *
from modules.ProgressBar  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''




def mutationCall(sortedLocations,cov,diffRate,prefix,species,VERBOSE):
  
   
    
    progressBar = ProgressBar(sys.stdout,"Calling Mutations...",'.','Complete',1,VERBOSE)
    progressBar.increment()


    mutations = MutationRecord(sortedLocations,prefix,cov,diffRate,species)

    while mutations.open:
       
        mutations.geneCandSearch()
        mutations.geneCandCall()
    
    progressBar.complete()

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-c", "--coverage", default = 5,  type=int, help="Minimum Coverage for Mutation")
    parser.add_option("-d", "--diffRate", default = 0.4, type=float, help="Minimum difference rate for mutation")
    parser.add_option("-s", "--chooseSpecies", default = "HUMAN", type='string', help="Change Species")
    parser.add_option("-p", "--prefix", default = "test", type='string', help="prefix file name")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    
    (options, args) = parser.parse_args()

    if options.chooseSpecies.upper() in ["MONKEY","RHESUS"]:          SPECIES="RHESUS"
    else:                                                           SPECIES="HUMAN"



    if len(args) == 1:
        mutationCall(args[0],options.coverage,options.diffRate,options.prefix,options.chooseSpecies,options.verbose)
    else:

        parser.print_help()
        print ""
        print "Example Usage: ./call_mutatations gene_map.srt --coverage 10 --diffRate 0.3 --verbose"
        sys.exit(2)

