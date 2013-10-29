#!/usr/bin/env python


import sys
import os
from collections import defaultdict as dd	
from modules.Utilities import *
from modules.SpliceRecord import *
from modules.GtfFile import *
from math import fabs




def fillUpExpression(myfile):
    myDict=dd(lambda: [0,0])
    if not myfile:
        return myDict
    for line in open(myfile):
        line=line.split()
        if len(line)!=3: continue
        myDict[line[0]][0] += int(line[1])
        myDict[line[0]][1] += int(line[2])
        splitName=line[0].split(",")
        if len(splitName) > 1:
            for n in splitName:
                myDict[n][1] += int(line[1])
    return myDict


def fillUpSplice(myfile):
    myDict = dd(lambda: [0,[]])
    if not myfile: return myDict
    for line in open(myfile):
        line=line.split()
        feature = line[0].split("|")
        splice = (feature[0],feature[-1]); data=tuple(feature[0:len(feature)-1])
        myDict[splice][0]+=int(line[1])
        myDict[splice][1].append(data)
        
    return myDict




def process_files(INSTRUCTIONS,exonFile,intronFile,spliceFile,genomeFile):

    if INSTRUCTIONS=="EXPRESSION":
        exTable=fillUpExpression(exonFile)
        intTable=fillUpExpression(intronFile)
        spliceTable=fillUpExpression(spliceFile)
        ALL_GENES = set( exTable.keys() + intTable.keys() + spliceTable.keys() )
        for s in ALL_GENES: 
            cnts = [exTable[s][0],exTable[s][1],intTable[s][0],intTable[s][1],spliceTable[s][0],spliceTable[s][1]]
            print s,"UNIQ/AMBIG-CNTS: EXONIC",cnts[0],cnts[1],"INTRONIC",cnts[2],cnts[3],"NOVEL-SPLICED",cnts[4],cnts[5] 

    elif INSTRUCTIONS == "SPLICEDATA":
        
        exTable=fillUpSplice(exonFile)
        intTable=fillUpSplice(intronFile)
        spliceTable=fillUpSplice(spliceFile)
    
        ALL_SITES = set( exTable.keys() + intTable.keys() + spliceTable.keys() )

        for s in ALL_SITES:
            cnts = exTable[s][0]+intTable[s][0]+spliceTable[s][0]
            data = exTable[s][1]+intTable[s][1]+spliceTable[s][1]
            #print data
            if len(data) == 1 or len(set(data))==1:
                print "|".join([p for p in data[0]])+"|"+s[1],cnts
            else:
                #print data
                CANONS=[d for d in data if d[4] == "exCANON_CANON_exCANON"]
                if len(CANONS)>0:
                    print "|".join([p for p in CANONS[0]])+"|"+s[1],cnts
                else:
                    data.sort(key= lambda x: len(x[1]),reverse=True)
                    print "|".join([p for p in data[0]])+"|"+s[1],cnts
    
    elif INSTRUCTIONS == "STATS":

        for f in [exonFile,intronFile,spliceFile,genomeFile]:
            for line in open(f):
                print line.strip()
            print ""

    else:
        errorQuit("INVALID INSTRUCTIONS")


if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    #parser.add_option("-k", "--key", action = 'store', default = None, help="phenotype key")
    parser.add_option("-e", "--exonFile", default = None, type='string', help="Path to exon mapping file")
    parser.add_option("-i", "--intronFile", default = None, type='string', help="Path to intron file")
    parser.add_option("-s", "--spliceFile", default = None, type='string', help="Path to splice file")
    parser.add_option("-g", "--genomeFile", default = None, type='string', help="Path to genome file")


    (options, args) = parser.parse_args()

    if len(args)!=1:
        print "NEED A ARGUMENT" 
        print args
        sys.exit()
   
    #elif not options.key:
    #    print "A key is required"
    #    sys.exit()
    else:
        process_files(args[0],options.exonFile,options.intronFile,options.spliceFile,options.genomeFile)
