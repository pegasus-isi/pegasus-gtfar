#!/usr/bin/env python


import sys
from modules.sam.SamFile import *
import os
				
def process_file(sam_mapping,keyFile,prefix,EXON,INTRON,strandSpecific):
    

    sam = SamFile(sam_mapping,prefix,strandSpecific)
    sam.addKey(keyFile,EXON,INTRON)

    while sam.open:
        sam.getRead()
        sam.processRead()
        sam.writeToBam()
        sam.writeGeneLocations()
    rpkmOut=open(prefix+'.rpkm','w');  spliceOut=open(prefix+'.spliceSites','w');  featOut=open(prefix+'.features','w')
    sam.printRPKMs(rpkmOut);           sam.printSplices(spliceOut);                sam.printFeatures(featOut)
    rpkmOut.close();                   spliceOut.close();                          featOut.close()
    sam.readOut.close();               sam.samOut.close();
    
    systemCall="sort -k11n,11 -k6,6 < "+prefix+"_gene.locs > "+prefix+"_gene.srtLocs"
    os.system(systemCall)
    print sam.index,'read mappings analyzed'
        
    
            











if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", "--key", default = None, type='string', help="Path to key file")
    parser.add_option("-p", "--prefix", default = 'foobar', type='string', help="Prefix for OutPut")
    parser.add_option("-e", "--exonic", action = 'store_true', default = False, help="Exonic only  input")
    parser.add_option("-i", "--intronic", action = 'store_true', default = False, help="Intronic Only input")
    parser.add_option("-s", "--strand", default = None, type='string', help="0,+, OR 16,-")

    (options, args) = parser.parse_args()

    if len(args)==0:
        print "PARSER FOR SAM FILE (FRESH FILE)"
        print ""
        print "USAGE: ./parse_sam_mapping.py MAP_FILE.sam -k KEYFILE.txt/pickle -p PREFIX (for output)"
    elif options.exonic== False and options.intronic == False:
        print "Intronic or Exonic Specification Required"
        sys.exit()
    elif len(args)==1:
        process_file(args[0],options.key,options.prefix,options.exonic,options.intronic,'16')
    else:
        print "TOO MANY ARGS" 
