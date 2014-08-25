#!/usr/bin/env python


import os
import sys
#from modules.fastq.FastqFile import *
import modules.ToolSet as gtTools 
#from modules.FastqFilter import FastqReads
from modules.SamFile import SamFile
				
                
if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Read Length")
    parser.add_option("-t", "--trimLengths", default = '50,75', type='string', help="Trimmable Lengths")
    parser.add_option("-q", "--lowqual", default = 1, type='int', help="Low Quality Score")
    parser.add_option("-a", "--avgqual", default = 15, type='int', help="Minimum Avg Quality Score")
    parser.add_option("-p", "--prefix", default = 'foo', type='string', help="Output Filename Prefix")
    parser.add_option("-s", "--strictness", default = 1, type='int', help="Output Filename Prefix")
    #parser.add_option("-m", "--multiple", default = 25, type='int', help="Trimming Multiple")

    (options, args) = parser.parse_args()

    try:                FILE=gtTools.fileGrab(args,0)
    except IndexError:
        parser.print_help()
        sys.exit()

    data  =  SamFile(FILE,options.prefix)
    while data.fileOpen:
        data.getSamReadData()
    data.printResults()
