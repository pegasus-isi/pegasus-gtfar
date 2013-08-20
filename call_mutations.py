#!/usr/bin/env python

import sys
#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/PIPELINE_v2/gtFar/python_src/gtfar_modules')
#from modules.gtf.GtfFile  import *
#from modules.mutations.MutationRecord  import *

from modules.GtfFile  import *
from modules.MutationRecord  import *
from modules.ProgressBar  import *
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



def mutationCall(sortedLocations,gtfFile,FTYPE,genomePath,readlen,cov,diffRate,prefix,ANNOTATE,VERBOSE):
  
   
    
    progressBar = ProgressBar(sys.stdout,"Calling Mutations...",'.','Complete',1,VERBOSE)
    progressBar.increment()


    mutations = MutationRecord(sortedLocations,prefix,FTYPE,cov,diffRate)
    ## Initialize a mutation record; pointer to a sorted file; candidate data struct; mutation sites data structs ##
    
    mutations.findCands()
    ## Iterate through mapping file to count substitutions for each position ##

    gtf = GtfFile(gtfFile,prefix,readlen,FTYPE,False) 
    ## Initialize a gtfFile class; class which holds all information from gtf file ##


    while gtf.open:

        progressBar.increment()

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
            gtf.printAnnotation(FTYPE)
        gtf.startNextChromosome()
    
    progressBar.complete()

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-l", "--readlength", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genomePath", default = None, type='string', help="Path to genome chr fasta files")
    parser.add_option("-c", "--coverage", default = 5,  type=int, help="Minimum Coverage for Mutation")
    parser.add_option("-d", "--diffRate", default = 0.5, type=float, help="Minimum difference rate for mutation")
    parser.add_option("-a", "--annotate", action = 'store_true', default = False,help="reannoate reference")
    
    parser.add_option("-r", "--ref", default = "NA", type='string', help="Reference Type Used for Mapping") 
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    
    (options, args) = parser.parse_args()

    if options.ref.upper()   in ["TRANSCRIPTS","EXONS","EXONIC"]:   REFTYPE="EXONIC"
    elif options.ref.upper() in ["INTRONS","INTRONIC"]:             REFTYPE="INTRONIC"
    elif options.ref.upper() in ["GENOME","GENOMIC","HG19"]:        REFTYPE="HG19"
    elif options.ref.upper() in ["GENES"]:                          REFTYPE="GENES"
    else:                                                           REFTYPE=None


    if len(args) == 2 and options.genomePath != None and REFTYPE != None and options.prefix != None:
        mutationCall(args[0],args[1],REFTYPE,options.genomePath,options.readlength,options.coverage,options.diffRate,options.prefix,options.annotate,options.verbose)
    else:

        parser.print_help()
        print ""
        print "Example Usage: ./call_mutatations gene_map.srt gencode.gtf --ref=exonic -g GENOME_PATH -l 100 --coverage 10 --diffRate 0.3 --prefix exonic --annotate --verbose"
        sys.exit(2)

