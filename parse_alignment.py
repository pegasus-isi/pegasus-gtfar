#!/usr/bin/env python


import sys
#from modules.MapFile2 import *
from modules.MapLine import *
from modules.MapData import *
from modules.ProgressBar import *
import os
				
def process_file(mapFile,keyFile,prefix,strandSpecific,VERBOSE):

    progressBar = ProgressBar(sys.stderr,"Parsing Alignment...",'.','Complete',100000,VERBOSE)

    mapLines = MapLines(mapFile)

    mapReads = MapRead(mapLines,keyFile,strandSpecific)
  
    mapData  = MapData(prefix,keyFile,mapLines.refType,mapLines.fileName)


    while mapLines.open:
        mapReads.loadRead()
        mapData.process(mapReads)
    mapData.printResult()


    progressBar.complete()


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", "--key", dest='key',default = None, type='string', help="    Path to key file")
    parser.add_option("-p", "--prefix", dest='prefix', default = None, type='string', help="   Prefix for OutPut")
    #parser.add_option("-r", "--ref", default = "NA", type='string', help="Reference Type Used for Mapping") 
    parser.add_option("-s", "--strandRule", default = "NA", type='string', help="Opposite,Match")
    #parser.add_option("-c", "--chooseSpecies", default = "HUMAN", type='string', help="Change Species")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    parser.add_option("-j", "--justExpression", action = 'store_true', default = False,  help="verbose output")



    (options, args) = parser.parse_args()


#    if options.ref.upper()   in ["TRANSCRIPTS","EXONS","EXONIC"]:   REFTYPE="EXONIC"
#    elif options.ref.upper() in ["INTRONS","INTRONIC"]:             REFTYPE="INTRONIC"
#    elif options.ref.upper() in ["GENOME","GENOMIC","HG19","HG"]:        REFTYPE="HG19"
#    elif options.ref.upper() in ["INTERGENIC"]:                     REFTYPE="INTERGENIC"
#    elif options.ref.upper() in ["GENES"]:                          REFTYPE="GENES"
#    else:                                                           REFTYPE=None

#    if options.chooseSpecies.upper() in ["MONKEY","RHESUS"]:          SPECIES="RHESUS"
#    else:                                                           SPECIES="HUMAN"

    if options.strandRule.upper() in ["MATCH","SAME"]:                          STRANDMATCH = "MATCH"
    elif options.strandRule.upper() in ["REVERSE","DIFFERENT","OPPOSITE"]:      STRANDMATCH = "OPPOSITE"
    else:                                                                       STRANDMATCH = None



    if len(args)==1  and options.key != None and options.prefix != None:
        process_file(args[0],options.key,options.prefix,STRANDMATCH,options.verbose)

    else:
        parser.print_help()
        print ""
        print "Example Usage: ./parse_alignment.py mymapping.sam -k gtf.key -p myexonicoutput --strand match"
        sys.exit(2)

    
 
