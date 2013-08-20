#!/usr/bin/env python


import sys
#from modules.file_types.MapFile import *
from modules.MapFile import *
from modules.ProgressBar import *
import os
				
def process_file(mapFile,keyFile,FTYPE,GAPPED,prefix,strandSpecific,VERBOSE):

    progressBar = ProgressBar(sys.stdout,"Parsing Alignment...",'.','Complete',1000000,VERBOSE)

    mapping = MapFile(mapFile,keyFile,FTYPE,GAPPED,prefix,strandSpecific)

    if mapping.invalid != False: return mapping.invalid
    

    while mapping.open:
        progressBar.increment()
        mapping.getReads()

        if FTYPE != "HG19":
            mapping.storeExpression()
            mapping.writeLocations()

    
    if FTYPE != "HG19":
        mapping.writeExpression()
        mapping.close()
        systemCall="sort -k14n,14 -k6,6 -k8n,8 < "+prefix+"_gene.loc > "+prefix+"_gene.srt"
        os.system(systemCall)
    else:
        mapping.writeNovelGenes()
    progressBar.complete()
    return "PASS"


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", "--key", dest='key',default = None, type='string', help="    Path to key file")
    parser.add_option("-p", "--prefix", dest='prefix', default = None, type='string', help="   Prefix for OutPut")
    parser.add_option("-r", "--ref", default = "NA", type='string', help="Reference Type Used for Mapping") 
    parser.add_option("-g", "--gapped", action = 'store_true', default = False, help="Flag to denote gapped alignment")
    parser.add_option("-s", "--strand", default = "NA", type='string', help="0,+, OR 16,-")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")



    (options, args) = parser.parse_args()

    PROGRAM_RESULT="PASS"

    if options.ref.upper()   in ["TRANSCRIPTS","EXONS","EXONIC"]:   REFTYPE="EXONIC"
    elif options.ref.upper() in ["INTRONS","INTRONIC"]:             REFTYPE="INTRONIC"
    elif options.ref.upper() in ["GENOME","GENOMIC","HG19","HG"]:        REFTYPE="HG19"
    elif options.ref.upper() in ["GENES"]:                          REFTYPE="GENES"
    else:                                                           REFTYPE=None



    if len(args)==1 and REFTYPE != None and options.key != None and options.prefix != None:
        if options.strand.upper() in [ "BOTH", "NA", "NONE", "EITHER" ]: options.strand = None
        PROGRAM_RESULT=process_file(args[0],options.key,REFTYPE,options.gapped,options.prefix,options.strand,options.verbose)

    else:
        parser.print_help()
        print ""
        print "Example Usage: ./parse_alignment.py mymapping.sam -r exonic -k gtf.key -p myexonicoutput --strand + --gapped"
        sys.exit(2)

    if PROGRAM_RESULT != "PASS":


        if PROGRAM_RESULT == "EMPTY_MAPFILE":
            print "WARNING: Empty Alignment File"
            sys.exit(0)
        else:
            if PROGRAM_RESULT == "MAPFILE_EXTENSION":
                print "ERROR: Mapping or Sam Alignment required (extension must be .mapping or .sam)\n"
            parser.print_help()
            print "Example Usage: ./parse_alignment.py mymapping.sam -r exonic -k gtf.key -p myexonicoutput --strand + --gapped"
            sys.exit(2)
    
        

