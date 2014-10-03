#!/usr/bin/env python


import sys
from modules.ClipRead import *
import os
				

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose", action = 'store_true', default = False,  help="verbose output")
    parser.add_option("-t", "--tag", default = 'gtfar', type='string',  help="tag data")



    (options, args) = parser.parse_args()




    if len(args)==1:
        try:
            fileType = args[0].split(".")[-1]
            fileHandle = open(args[0])
            ClipReads = ClipReads(fileHandle,fileType,options.tag)
            while ClipReads.fileOpen:
                ClipReads.getNextRead()
            ClipReads.printData()
        except IOError:
            sys.exit() 

    else:
        parser.print_help()
        print ""
        print "Example Usage: ./parse_alignment.py mymapping.sam -p myexonicoutput"
        sys.exit(2)

    
 
