#!/usr/bin/env python


import sys
from modules.file_types.MapFile import *
import os
				
def process_file(mapFile,keyFile,prefix,FTYPE,strandSpecific):
   
    mapping = MapFile(mapFile,prefix,FTYPE,strandSpecific)

    mapping.loadKey(keyFile)
    
    while mapping.open:
        mapping.getReads()
        if FTYPE != "HG19":
            mapping.storeExpression()
            mapping.writeLocations()
    
    if FTYPE != "HG19":
        mapping.writeExpression()
        #print mapping.index,'read maps analyzed'
        mapping.close()
        systemCall="sort -k14n,14 -k6,6 -k8n,8 < "+prefix+"_gene.loc > "+prefix+"_gene.srt"
        os.system(systemCall)
        sys.exit()
    else:
        mapping.writeNovelGenes()


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", "--key", default = None, type='string', help="Path to key file")
    parser.add_option("-p", "--prefix", default = 'foobar', type='string', help="Prefix for OutPut")
    parser.add_option("-e", "--exonic", action = 'store_true', default = False, help="Exonic only  input")
    parser.add_option("-i", "--intronic", action = 'store_true', default = False, help="Intronic Only input")
    parser.add_option("-g", "--gapped", action = 'store_true', default = False, help="Gapped input")
    parser.add_option("-z", "--hg19", action = 'store_true', default = False, help="Gapped input")
    parser.add_option("-s", "--strand", default = None, type='string', help="0,+, OR 16,-")

    (options, args) = parser.parse_args()

    if len(args)==0:
        print "PARSER FOR ALIGNMENT FILE (FRESH FILE)"
        print ""
        print "USAGE: ./parse_alignment_mapping.py file.mapping --exonic -k KEYFILE.txt/pickle -p PREFIX (for output)"
    elif len(args)==1:
        mapFile=args[0]

        if options.gapped == True and options.hg19 == False:
            process_file(args[0],options.key,options.prefix,"GAPPED",'+')
        elif options.exonic == True and options.intronic == False:
            process_file(args[0],options.key,options.prefix,"EXONIC",'+')
        elif options.exonic == False and options.intronic == True:
            process_file(args[0],options.key,options.prefix,"INTRONIC",'+')
        elif options.exonic == False and options.intronic == False and options.hg19 == True:
            process_file(args[0],options.key,options.prefix,"HG19",'+')
        else:
            print "Intronic/Gapped/Exonic Specification Required"
            sys.exit()
    else:
        print "TOO MANY ARGS"
        sys.exit()
