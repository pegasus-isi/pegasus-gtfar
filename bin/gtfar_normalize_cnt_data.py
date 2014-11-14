#!/usr/bin/env python

import sys
import os
from collections import defaultdict as dd
'''
This program requires a gencode annotation file (gtf) and a path to the chr fasta files referenced in the gtf file
example: HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19
'''
mito = dd(bool)
mito.update({'NR_003285.2': True,  'NR_003286.1': True, 'NR_003287.1' : True, 'NR_MITO' : True})
#mito = dict(dd(bool).items() + {'NR_003285.2': True,  'NR_003286.1': True, 'NR_003287.1' : True, 'NR_MITO' : True}.items())

gene_cnts = dd(lambda: dd(int))
total_cnts = dd(int)
total_ranks = dd(int)

def normalize(files,prefix,MITO=False,AMBIG=False,RPM=True,QUANTILE=True,RPM_factor=10000000.0):
    
    if not AMBIG:
        for f in files:
            total = 0
            for line in open(f):
                line = line.split()
                name,val=line[0],int(line[1])
                if mito[name.split(',')[0]]:
                    if not MITO: continue
                gene_cnts[name][f] = val
                total+=val
            total_cnts[f] = total
            
    genes = sorted(gene_cnts.keys())
    
    if RPM:
        MEAN_TOTAL = sum(total_cnts.values()) / float(len(total_cnts.keys()))
        t_len = len(str(int(MEAN_TOTAL)))-1
        RPM_FACTOR = round(MEAN_TOTAL,-1*t_len)
        
        genes = sorted(gene_cnts.keys())

        for filename in total_cnts.keys():
            if prefix == None:
                outname = ".".join([x for x in filename.split(".")[0:len(filename.split("."))-1]])+".rtn"
                out=open(outname,"w")
            factor = RPM_FACTOR / total_cnts[filename]
            for g in genes:
                out.write("%s %s\n" % (g,gene_cnts[g][filename]*factor))
    if QUANTILE:
        file_num = float(len(total_cnts.keys()))
        for filename in total_cnts.keys():
            value_sort = sorted([ (gene_cnts[g][filename],g) for g in genes ])
            for i in range(len(value_sort)):
                total_ranks[i]+=value_sort[i][0]
                gene_cnts[value_sort[i][1]][filename] = i
        for filename in total_cnts.keys():
            if prefix == None:
                outname = ".".join([x for x in filename.split(".")[0:len(filename.split("."))-1]])+".quantile"
                out=open(outname,"w")
            for g in genes:
                out.write("%s %s\n" % (g,total_ranks[gene_cnts[g][filename]]/file_num ))




if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-l", "--includeMito",action='store_true', help="Ignore Mito/Ribo Cnts")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")
    parser.add_option("--quantile", action='store_true', help="Quantile Normalization")
    parser.add_option("--rpm", action='store_true', help="RPM normalization")
    parser.add_option("--ambiguous", action='store_true', help="Include Ambiguous Reads")

    #parser.add_option("-g", "--genome", default = None, type='string', help="Path to genome autosome fasta files")
    #parser.add_option("-c", "--config", default = None, type='string', help="Path to configuration file")
    #parser.add_option("-k", "--key", default = None, type='string', help="Path to key file")
    #parser.add_option("-f", "--filters", default = None, type='string', help="Path to filter file")

    (options, args) = parser.parse_args()



if not options.quantile and not options.rpm:
    print "NEED AT LEAST ONE NORMALIZATION TYPE (--quantile or --rpm)"
    sys.exit()
elif len(args) == 0:
    print "NEED AT LEAST ONE FILE TO NORMALIZE"
    sys.exit()
else:
    normalize(args,options.prefix,options.includeMito,options.ambiguous,options.rpm,options.quantile)



