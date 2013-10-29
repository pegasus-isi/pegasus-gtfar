#!/usr/bin/env python


import sys
from random import randrange
from random import choice
from collections import defaultdict as dd
from modules.Utilities import *
from math import fabs

##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class SpliceRecord:
    def __init__(self):
        
        self.geneJumps = dd(lambda: dd(int))
            
        self.validJxns = dd(list)


    def addKey(self, key):
        self.geneKey={}; self.geneExonTable = dd(list); self.geneJxnTable =dd(list); self.geneSplitTable = dd(lambda: dd(list))
        for line in open(key):
            line=line.split(); feature=line[0].split("|");  fName=feature[0]; fType=feature[len(feature)-1]
            if fType == "GENE":
                self.geneKey[fName]=line[0]
                self.geneExonTable[fName]=[(int(x.split("-")[0]),int(x.split("-")[1])) for x in line[1].split(",")]
            elif fType == "JXNS":
                if len(line)>1:
                    self.geneJxnTable[fName] = [(int(x.split(",")[0]),int(x.split(",")[1])) for x in line[1].split("|")]
                    tmpSplits = sorted([(g[1],g) for g in self.geneJxnTable[fName]]+[(g[0],g) for g in self.geneJxnTable[fName]])
                    tmpExons=sorted(self.geneExonTable[fName])
                    k=0
                    for i in range(len(tmpSplits)):
                        while tmpSplits[i][0]>tmpExons[k][1]:   k+=1
                        if tmpSplits[i][0] <= tmpExons[k][1] and tmpSplits[i][0] >= tmpExons[k][0]:
                           self.geneSplitTable[fName][tmpExons[k]].append(tmpSplits[i])
                        else:
                            print "FUCK", i,tmpSplits[i],tmpExons[k-1],tmpExons[k]
                            sys.exit()


    @staticmethod
    def findNearestSplice(site,type,key):
        if type == "END":
            ENDS = [x for x in key if x[0] == x[1][0]]
            if len(ENDS)==0: return ["NOVEL",site,None,None]
            if site in [x[0] for x in ENDS]:
                cands = [x[1][1] for x in ENDS if x[0] == site]
                return ["CANONICAL",site,"PREDICTED-PARTNERS:"]+cands
            else:
                neighbors = sorted([(abs(x[0]-site),x) for x in ENDS])
                dist = neighbors[0][0]; spot=neighbors[0][1][0]
                cands = [x[1][1] for x in ENDS if x[0] == spot]
                return ["NEIGHBOR",site,spot,dist,"NB-PARTNERS:"]+cands
        else:
            STARTS = [ x for x in key if x[0] == x[1][1]]
            if len(STARTS) == 0: return ["NOVEL",site,None,None]
            if site in [x[0] for x in STARTS]:
                cands = [x[1][0] for x in STARTS if x[0] == site]
                return ["CANONICAL",site,"PREDICTED-PARTNERS:"]+cands
            else:
                neighbors = sorted([(abs(x[0]-site),x) for x in STARTS])
                dist = neighbors[0][0]; spot=neighbors[0][1][0]
                cands = [x[1][0] for x in STARTS if x[0] == spot]
                return ["NEIGHBOR",site,spot,dist,"NB-PARTNERS:"]+cands

    @staticmethod
    def categorizeSplices(endData,startData,cnt,geneStarts,geneEnds):
        JD = startData[1]-endData[1]
        myCands = [[endData[1],startData[1],cnt]]

        if endData[0] == "NEIGHBOR" and startData[0] == "NEIGHBOR":
            if endData[2] in startData[5::] and ((endData[1]-endData[2]) == (startData[1]-startData[2])): myCands.append([endData[2],startData[2],cnt])
            elif endData[3] < 100 and endData[2]+JD in geneStarts and endData[2]+JD not in endData[5::]:       myCands.append([endData[2],endData[2]+JD,cnt])
            elif startData[3] < 100 and startData[2]-JD in geneEnds and startData[2]-JD not in startData[5::]: myCands.append([startData[2]-JD,startData[2],cnt])

        if endData[0] == "NEIGHBOR" and startData[0] == "NOVEL":
            if endData[2] < endData[1] and endData[2] < (startData[1] - endData[3]): myCands.append([endData[2],startData[1]-endData[3],cnt])
            elif endData[2] > endData[1] and endData[2] < (startData[1] + endData[3]):  myCands.append([endData[2],startData[1]+endData[3],cnt])

        elif endData[0] == "NOVEL" and startData[0] == "NEIGHBOR":
            if startData[2] < startData[1] and startData[2] > (endData[1] - startData[3]):      myCands.append([endData[1]-startData[3],startData[2],cnt])
            elif startData[2] > startData[1] and startData[2] > (endData[1] + startData[3]):    myCands.append([endData[1]+startData[3],startData[2],cnt])
            
        return myCands 





    def addGappedAlignmentFiles(self, args):
        gene_jumps = dd(lambda: dd(int))
        for f in args:
            fname=f.split(".")[0]
            for line in open(f):
                line=line.split()
                feat=line[2]; spot=int(line[3]); cigar=line[5]
                gene_jumps[feat.split("|")[0]][tuple(cigarToLoc(spot,cigar)[1:3])]+=1
        
        for gene in gene_jumps:
            geneExons=self.geneExonTable[gene];         geneJxns=self.geneJxnTable[gene];       geneSplits = self.geneSplitTable[gene]
            geneStarts=set([x[1] for x in geneJxns]);   geneEnds=set([x[0] for x in geneJxns])
        
            jxns = sorted([((x[0]+(geneExons[0][0]-1),x[1]+(geneExons[0][0]-1)),gene_jumps[gene][x]) for x in gene_jumps[gene]]); k=0; spliceCands = []
    
            for j in jxns:
                endInfo = []; startInfo = []; jEnd =j[0][0]; jStart =j[0][1]; jumpDist = jStart-jEnd; JD=jumpDist
                while k < len(geneExons)-1 and jEnd > geneExons[k][1]: k+=1
                m=k
                if jEnd < geneExons[k][0]:                                  endInfo = ["INTRON",geneExons[k-1],jEnd-geneExons[k-1][1]]
                elif geneExons[k][0] <= jEnd and jEnd <= geneExons[k][1]:   endInfo = ["EXON",geneExons[k]]
                else:                                                       errorQuit("WRID")
            
                while m<len(geneExons) and jStart > geneExons[m][1]: m+=1
            
                if jStart < geneExons[m][0]:                                     startInfo = ["INTRON",geneExons[m],geneExons[m][0]-jStart]
                elif geneExons[m][0] <= jStart and jStart <= geneExons[m][1]:    startInfo = ["EXON",geneExons[m]]
                else:                                                            errorQuit("WigRID")
            
                endData,startData   =   self.findNearestSplice(jEnd,"END",geneSplits[endInfo[1]]), self.findNearestSplice(jStart,"START",geneSplits[startInfo[1]])
                spliceCands += self.categorizeSplices(endData,startData,j[1],geneStarts,geneEnds)

            spliceCands.sort()
            splicePass = []; lastCand = spliceCands[0]; lastSpot = lastCand[0:2]; lastCnt=lastCand[2]
            for s in spliceCands[1::]:
                if s[0:2] == lastSpot: lastCnt+=s[2]
                elif lastCnt>1: splicePass.append(lastSpot)
                lastCand=s; lastSpot=lastCand[0:2]; lastCnt=lastCand[2]
            if lastCnt>1: splicePass.append(lastSpot)

            self.validJxns[gene] = sorted(splicePass)
                            
