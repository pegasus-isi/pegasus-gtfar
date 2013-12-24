#!/usr/bin/env python


import sys
import os
from collections import defaultdict as dd	
import modules.SampleAnalysis #import *
from modules.Utilities import *

def median(mylist):
    return mylist[len(mylist)/2]

def quartile(mylist):
    return mylist[3*(len(mylist) /4) ]


def pearsonProcess(args):
    GENE_CNTS=dd(lambda: [0,0])
    k=0
    sampleData = modules.SampleAnalysis.SampleData()
    for f in args:
        sampleData.addMember(modules.SampleAnalysis.Sample(f,"LOAD_GENE_CNTS"))
        #samples.addMember(modules.SampleAnalysis.Sample(f,"LOAD_GENE_CNTS"))
   
    sampleData.mergeGenes()

    sampleData.pairwiseCorrelation(0,1)


    

        













def process_files(args,EXPR,key,NORM):

    outcome = {}; samples = []; genes = dd(lambda: 1.0)
    for line in open(key):
        line=line.split()
        outcome[line[0].split('.')[0]] = line[1]
    
    if EXPR:
        geneTable = dd(lambda: dd(float)); 
        exNormal, intNormal, cntTable = {},{},{}
        argCnt =0; fNames = []; multiGenes = dd(list)

        for f in args:
            fname=f.split(".")[0]
            mySample = SampleExpr(fname)
            for line in open(f):
                line=line.split()
                geneName= line[0]; tmpCnts = float(line[1])
                if tmpCnts > 0: 
                    mySample.addGeneCnt(geneName,tmpCnts)
                    genes[geneName]*=tmpCnts
                if len(geneName.split("/")) > 1:
                    gList = geneName.split("/")
                    for g in gList: multiGenes[g].append(geneName)
            samples.append(mySample)
        
        for s in samples:
            s.process(genes,len(samples))
            s.addResult(outcome[s.name])
            s.normalize(NORM)
        for g in genes:

            geneData = [(s.phenotype,s.norm[g],s.name) for s in samples]
            #intData = [(s.phenotype,s.norm[g][1],s.name) for s in samples]
            #print g,exData,intData
            myGene = GeneExpr(g,geneData)
            myGene.calculate_continuous_correlation()   
            
            
            

    else:
        spliceDict=dd(lambda: [None,0])
        for f in args:
            for line in open(f):
                line=line.split()
                if line[1] == "INTRONIC": continue
                elif line[1] == "EXONIC" or line[1]=="KNOWN-JXNS":
                    spliceDict[line[0]][0]="KNOWN"
                    spliceDict[line[0]][1]+= int(line[2])
                else:
                    if spliceDict[line[0]][0]:
                        spliceDict[line[0]][1] += int(line[2])
                    else:
                        spliceDict[line[0]][0] = "NOVEL"
                        spliceDict[line[0]][1] += int(line[2])
        for s in spliceDict:
            cnts=spliceDict[s]
            print s,cnts[0],cnts[1]


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("--twoFilePearson",  action = 'store_true', default = False, help="two file correlation")
    parser.add_option("-e", "--expr",  action = 'store_true', default = False, help="expression cnts")
    parser.add_option("-s", "--splice", action = 'store_true', default = False, help="splice cnts")
    parser.add_option("-k", "--key", action = 'store', default = None, help="phenotype key")
    parser.add_option("-n", "--norm", action = 'store', default = "MEAN", help="normalize")

    (options, args) = parser.parse_args()

    if options.twoFilePearson:
        if len(args) != 2:  errorQuit("Need two files for pearson check")
        pearsonProcess(args)
        sys.exit()

    if len(args)<2:
        print "WE NEED MORE FILES" 
        print args
        sys.exit()
    elif options.splice == options.expr:
        print "option expression or option splice required"
        sys.exit()
    elif not options.key:
        print "A key is required"
        sys.exit()
    else:
        NORM=options.norm.upper()
        process_files(args,options.expr,options.key,NORM)
