#!/usr/bin/env python


import sys
#from modules.MapFile2 import *
from modules.MapLine import *
from modules.MapData import *
from modules.ProgressBar import *
import os
				
def process_file(mapFile,prefix,strandSpecific,VERBOSE):

    progressBar = ProgressBar(sys.stderr,"Parsing Alignment...",'.','Complete',100000,VERBOSE)

    mapLines = MapLines(mapFile)


    mapReads = MapRead(mapLines,strandSpecific)  
    mapData  = MapData(prefix,mapLines.fileName)


    while mapLines.open:
        mapReads.loadRead()
        mapData.process(mapReads)
    mapData.printResult()


    progressBar.complete()


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--prefix", dest='prefix', default = None, type='string', help="   Prefix for OutPut")
    parser.add_option("-s", "--strandRule", default = "NA", type='string', help="Opposite,Match")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    parser.add_option("-j", "--justExpression", action = 'store_true', default = False,  help="verbose output")



    (options, args) = parser.parse_args()



    if options.strandRule.upper() in ["MATCH","SAME"]:                          STRANDMATCH = "MATCH"
    elif options.strandRule.upper() in ["REVERSE","DIFFERENT","OPPOSITE"]:      STRANDMATCH = "OPPOSITE"
    else:                                                                       STRANDMATCH = None



    if len(args)==1  and options.prefix != None:
        process_file(args[0],options.prefix,STRANDMATCH,options.verbose)

    else:
        parser.print_help()
        print ""
        print "Example Usage: ./parse_alignment.py mymapping.sam -p myexonicoutput --strand match"
        sys.exit(2)

    
 
