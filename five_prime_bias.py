#!/usr/bin/env python

import sys
#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/PIPELINE_v2/gtFar/python_src/gtfar_modules')
from modules.gtf.GtfFile  import *
from modules.mutations.MutationRecord  import *

'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



def mutationCall(sortedLocations,fName,genomePath,prefix,TYPE,readlen,coverage,diffRate,silent,printMuts,ANNOTATE):
   


    mutations = MutationRecord(sortedLocations,prefix,TYPE,coverage,diffRate)
    ## Initialize a mutation record; pointer to a sorted file; candidate data struct; mutation sites data structs ##
    
    mutations.findCands()
    ## Iterate through mapping file to count substitutions for each position ##

    gtf = GtfFile(fName,prefix,readlen,TYPE,False) 
    ## Initialize a gtfFile class; class which holds all information from gtf file ##


    while gtf.open:

        gtf.loadGenesOnChromosome();    gtf.uniquifySeqs(SILENT=True)
        ## Load the genes corresponding to the chromosome ## 
        
        while mutations.chr == gtf.chr:
            #print "Calling Chr Cands",mutations.chr,TYPE
            mutations.recordChrCands()
            ## Add detailed information for each read which alignes to a candidate ##
            #print "Evaluating Chr Cands",gtf.chr
            gtf.addFasta(genomePath+'/'+gtf.chr+'.fa')
            ## Load Fasta Sequence for the Chromosome ##
            gtf.evaluateMutations(mutations.chrSites)
            ## Using Fasta Sequence and Junction Information; whether or not a mutation exists ##
        if ANNOTATE:
            gtf.printAnnotation(TYPE)
        gtf.startNextChromosome()

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mutations', type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genomePath", default = None, type='string', help="Path to genome chr fasta files")
    parser.add_option("-c", "--coverage", default = 5,  type=int, help="Minimum Coverage for Mutation")
    parser.add_option("-d", "--diffRate", default = 0.5, type=float, help="Minimum difference rate for mutation")
    parser.add_option("-s", "--silent", action = 'store_true', default = False,help="silent the output")
    parser.add_option("-m", "--printMutations", action = 'store_true', default = False,help="print mutations")
    parser.add_option("-a", "--annotate", action = 'store_true', default = True,help="print mutations")
    
    
    
    parser.add_option("-e", "--exonic", action = 'store_true', default = False, help="Exonic only  input")
    parser.add_option("-i", "--intronic", action = 'store_true', default = False, help="Intronic only  input")
    (options, args) = parser.parse_args()


    if len(args)!=2:
        print "TWO ARGS REQUIRED"
        print "Usage: ./call_mutations_from_locations.py sorted_mapfile.sam gencode_file.gtf -g PATH_TO_GENOME -p OUTPUT_PREFIX --exonic/intronic"
    else:
        if options.genomePath == None:
            print "a genome path is required"
            print "Usage: ./call_mutations_from_locations.py sorted_mapfile.sam gencode_file.gtf -g PATH_TO_GENOME -p OUTPUT_PREFIX --exonic/intronic"
            sys.exit()
        if options.exonic:
            mutationCall(args[0],args[1],options.genomePath,options.prefix,'EXONIC',options.readlen,options.coverage,options.diffRate,options.silent,options.printMutations,options.annotate)
        elif options.intronic:
            mutationCall(args[0],args[1],options.genomePath,options.prefix,'INTRONIC',options.readlen,options.coverage,options.diffRate,options.silent,options.printMutations,options.annotate)
        else:

            print "NEED ONE"



