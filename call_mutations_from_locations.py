#!/usr/bin/env python

import sys
#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/PIPELINE_v2/gtFar/python_src/gtfar_modules')
from modules.gtf.GtfFile  import *
from modules.mutations.MutationRecord  import *

'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''



def mutationCall(sortedLocations,fName,genomePath,prefix,TYPE,readlen,coverage,diffRate):
   
    
    
    gtf = GtfFile(fName,prefix,readlen,printKEY=False)
    mutations = MutationRecord(sortedLocations,prefix,TYPE,coverage,diffRate)
    while gtf.open:
        gtf.loadGenesOnChromosome();  gtf.addFasta(genomePath+'/'+gtf.chr+'.fa');   gtf.uniquifySeqs(SILENT=True)
        print "OK STARTING MUTATIONS" 
        while mutations.chr == gtf.chr:
            print "CANDING"
            mutations.findCands(gtf);
            print "EVALIN"
            mutations.evalAndPrint();
            print "NEXT-GENE"
            mutations.nextGene()
            print "ANNOTATING"
        gtf.printAnnotation(TYPE)
        gtf.startNextChromosome()
    

    

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--readlen", default = 100, type='int', help="Expected Read Length")
    parser.add_option("-p", "--prefix", default = 'mutations', type='string', help="Output Filename Prefix")
    parser.add_option("-g", "--genomePath", default = None, type='string', help="Path to genome chr fasta files")
    parser.add_option("-c", "--coverage", default = 1, type=int, help="Minimum Coverage for Mutation")
    parser.add_option("-d", "--diffRate", default = 0.1, type=float, help="Minimum difference rate for mutation")
    
    
    
    parser.add_option("-e", "--exonic", action = 'store_true', default = False, help="Exonic only  input")
    parser.add_option("-i", "--intronic", action = 'store_true', default = False, help="Intronic only  input")
    (options, args) = parser.parse_args()


    if len(args)!=2:
        print "TWO ARGS REQUIRED"
        print "Usage: ./call_mutations_from_locations.py sorted_mapfile.sam gencode_file.gtf -g PATH_TO_GENOME -p OUTPUT_PREFIX"
    else:
        if options.exonic:
            mutationCall(args[0],args[1],options.genomePath,options.prefix,'EXONIC',options.readlen,options.coverage,options.diffRate)
        elif options.intronic:
            mutationCall(args[0],args[1],options.genomePath,options.prefix,'INTRONIC',options.readlen,options.coverage,options.diffRate)
        else:
            print "NEED ONE"



