#!/usr/bin/env python


import sys
#from MutationSite import *
from MutationCand import *
from Sequence import *
from Utilities import *


#from ..tools.gtTools import *

from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationRecord:
    def __init__(self,fileHandle,prefix,TYPE,coverage,diffRate):
        self.fileHandle = fileHandle
        self.open = False
        self.nextLine()
        self.fType = TYPE; self.prefix = prefix
        self.coverage_parameter = coverage;   self.rate_parameter = diffRate;     self.endLength=-1

###########################################################################################################################################################

    def nextLine(self):
        if self.open:
            self.line = self.fStream.readline().split()
            if len(self.line)>0:
                self.mapChr, self.mapGene = self.line[3], self.line[5]

                self.mapPos = [int(x) for x in self.line[7].split(",")]
                
                self.mapRead,self.mapRef,self.mapQual = self.line[9], self.line[10], self.line[11]
            else:
                self.mapGene = None
                self.mapChr  = None
                self.open = False
                self.fStream.close()
        else:
            self.fStream = open(self.fileHandle)
            self.open    = True
            self.line = self.fStream.readline().split()
            self.mapChr, self.mapGene = self.line[3], self.line[5]
            self.mapPos = [int(x) for x in self.line[7].split(",")]
            self.mapRead,self.mapRef,self.mapQual = self.line[9], self.line[10], self.line[11]
            ##########################
            self.chr  = self.mapChr
            self.gene = self.mapGene

###########################################################################################################################################################
            

    def findCands(self):
        self.cands = {}; self.candMut = 3; candCnt=0
        while self.open:
            geneCands = dd(int)
            while self.mapGene == self.gene:
            
                if int(self.line[12]) != 0 and self.line[1]=='UNIQ':

                    start = 0
                    spliceBuffer = 3
                    diffs=[]
                    for i in range(0,len(self.mapPos),2):
                        if self.mapPos[i+1]-self.mapPos[i] > spliceBuffer*2:
                            end = start + (self.mapPos[i+1] - self.mapPos[i] +1)
                            readSeq = self.mapRead[start:end];    refSeq  = self.mapRef[start:end]
                            diffs.extend([self.mapPos[i]+j for j in range(spliceBuffer,len(readSeq)-spliceBuffer) if readSeq[j]!=refSeq[j]])
                     #       if 1358 in diffs:
                     #           print self.line
                     #           print "YEAH"
                        start = start + (self.mapPos[i+1] - self.mapPos[i] + 1)
                    for d in diffs:
                        geneCands[d]+=1
                
                self.nextLine()
            
            passCands=[]
            for k in geneCands:
                if geneCands[k] >= self.candMut:
                    passCands.append(k); candCnt+=1
            self.cands[self.gene] = sorted(passCands,reverse=True)

            if not self.open:
                #print "Stored",candCnt,"Cands";
                self.nextLine()
                return 
            else:
                self.gene = self.mapGene

###########################################################################################################################################################
            
    def recordChrCands(self):

        self.chrSites = dd(list)
        ## ITERATE THROUGH SORTED FILE UNTIL THE LINES NO LONGER CORRESPOND TO THE CORRECT CHROMOSOME ##
        while self.mapChr == self.chr:

            ## ITERATE THROUGH FILE UNTIL A GENE WITH CANDIDATES IS FOUND ##
            while self.cands[self.mapGene] == []:
                self.nextLine()
                if not self.open:
                    break 
            ## LET FIRST OCCURANCE OF A CANDIDATE GENE BE MARKED ##
            if self.open:
                self.gene = self.mapGene
            
            ## IF THIS GENE CORRESPONDS TO A FUTURE CHROSOMOME; BREAK OUT OF LOOP #
                if self.mapChr != self.chr:
                    break
            ## INITIALIZE GENE CANDIDATES FOR THE GENE ##
                myCands = self.cands[self.gene]
                mySites = [MutationCand(x) for x in myCands]
            else:
                break

            ## AS LONG AS GENE CANDIDATES REMAIN - IF THE READ POSITION IS PAST THE CANDIDATE EVALUATE; OTHERWISE RECORD OR ITERATE THROUGH READS ##
            while len(myCands) > 0:

                if self.mapPos[0] > myCands[-1] or self.mapGene != self.gene:

                    ALTS = sum(mySites[-1].cnts) - max(mySites[-1].cnts)
                    TOT  = float(sum(mySites[-1].cnts))
                    


                    if TOT !=0 and (ALTS / TOT > 0.05 or mySites[-1].diffs / TOT > 0.05):
                        self.chrSites[self.gene].append(mySites[-1])
                    mySites.pop();  myCands.pop()

                else:
                    if self.mapPos[-1] >= myCands[-1] and self.line[1]!="AMBIG":

                        rType = self.line[1]; aType = self.line[2]
                        for i in range(len(myCands)-1,-1,-1):
                            start = 0
                            for j in range(0,len(self.mapPos),2):
                                if myCands[i] >= self.mapPos[j] and myCands[i] <= self.mapPos[j+1]:

                                    end = start+self.mapPos[j+1]-self.mapPos[j]+1;    candPos = myCands[i]-self.mapPos[j];  readIdx = start+candPos
                                    readBase = self.mapRead[start:end][candPos];  refBase = self.mapRef[start:end][candPos]; qualScr = self.mapQual[start:end][candPos]
                                    mySites[i].record(self.line[1],self.line[2],readIdx,readBase,refBase,qualScr)

                                start = start + (self.mapPos[j+1] - self.mapPos[j] + 1)

                            if self.mapPos[-1] < myCands[i]: break


                            #if not SPAN and self.mapPos[j+1] < 
                
                    self.nextLine()
            while self.mapGene == self.gene:
                self.nextLine()
        self.chr = self.mapChr
                
        
###########################################################################################################################################################
###########################################################################################################################################################
