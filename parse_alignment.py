#!/usr/bin/env python


import sys
from modules.MapRead import *
#from modules.ProgressBar import *
import os
				
def process_file(mapFile,strandSpecific,VERBOSE):

    #progressBar = ProgressBar(sys.stderr,"Parsing Alignment...",'.','Complete',100000,VERBOSE)

    mapReads = MapReads(mapFile,strandSpecific)

    while mapReads.fileOpen:
        mapReads.getNextRead()
        mapReads.printData()


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-s", "--strandRule", default = "NA", type='string', help="Opposite,Match")
    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    parser.add_option("-j", "--justExpression", action = 'store_true', default = False,  help="verbose output")



    (options, args) = parser.parse_args()



    if options.strandRule.upper() in ["MATCH","SAME"]:                          STRANDMATCH = "MATCH"
    elif options.strandRule.upper() in ["REVERSE","DIFFERENT","OPPOSITE"]:      STRANDMATCH = "OPPOSITE"
    else:                                                                       STRANDMATCH = None



    if len(args)==1:
        try:
            fileType = args[0].split(".")[-1]
            fileHandle = open(args[0])
            mapReads = MapReads(fileHandle,fileType,STRANDMATCH)
            while mapReads.fileOpen:
                mapReads.getNextRead()
                mapReads.printData()
        except IOError:
            sys.exit() 

    else:
        parser.print_help()
        print ""
        print "Example Usage: ./parse_alignment.py mymapping.sam -p myexonicoutput --strand match"
        sys.exit(2)

    
 
