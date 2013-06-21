#!/usr/bin/env python


import sys
from modules.sam.SamFile import *
				
def process_file(sam_mapping,keyFile,prefix):
    sam = SamFile(sam_mapping)
    sam.addKey(keyFile)
    samOut=open(prefix+'_viz.sam','w')
    readOut=open(prefix+'_gene.Locations','w')
    while sam.open:
        read=sam.getRead(strandSpecific='16')
        read.relocate()
        read.visualizeBAM(samOut,SHOWAMBIG=True)
        read.visualizeGeneLocation(readOut)
        sam.storeGeneExp(read)
    rpkmOut=open(prefix+'.rpkm','w')
    sam.printRPKMs(rpkmOut)
    print sam.index,'read mappings analyzed'
        
    
            











if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", "--key", default = None, type='string', help="Path to key file")
    parser.add_option("-p", "--prefix", default = 'foobar', type='string', help="Prefix for OutPut")

    (options, args) = parser.parse_args()

    process_file(args[0],options.key,options.prefix)

