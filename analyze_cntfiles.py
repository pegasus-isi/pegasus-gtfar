#!/usr/bin/env python


import os
import sys
#from modules.fastq.FastqFile import *
import modules.ToolSet as gtTools 
#from modules.FastqFilter import FastqReads
from modules.CntFiles import CntFiles
				
                
if __name__ == '__main__':

    from optparse import OptionParser
    from optparse import Option, OptionValueError
    usage = "usage: ./%prog [options] data_file"
    
    
    class MultipleOption(Option):
        ACTIONS = Option.ACTIONS + ("extend",)
        STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
        TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
        ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

        def take_action(self, action, dest, opt, value, values, parser):
            if action == "extend":
                values.ensure_value(dest, []).append(value)
            else:
                Option.take_action(self, action, dest, opt, value, values, parser)
    
    
    parser = OptionParser(usage=usage)
    #parser = OptionParser(option_class=MultipleOption,usage=usage)

    parser.add_option("-k", "--key", default = None, type='string', help="Key File")
    parser.add_option("-t", "--type", default = 'binary', type='string', help="binary or continuous")
    #parser.add_option("-q", "--lowqual", default = 1, type='int', help="Low Quality Score")
    #parser.add_option("-a", "--avgqual", default = 15, type='int', help="Minimum Avg Quality Score")
    parser.add_option("-p", "--prefix", default = 'foo', type='string', help="Output Filename Prefix")
    #parser.add_option("-s", "--strictness", default = 1, type='int', help="Output Filename Prefix")
    #parser.add_option("-m", "--multiple", default = 25, type='int', help="Trimming Multiple")

    (options, args) = parser.parse_args()

    #try:                FILE=gtTools.fileGrab(args,0)
    #except IndexError:
     #   parser.print_help()
      #  sys.exit()

    if options.key == None: 
        sys.stderr.write("No key supplied\n")
        parser.print_help()
        sys.exit()

    data = CntFiles(args,options.key,options.type,options.prefix)
    data.analyze()
    data.printResults()
    #data  =  SamFile(FILE,options.prefix)
    #while data.fileOpen:
    #    data.getSamReadData()
    #data.printResults()
