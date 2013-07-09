#!/usr/bin/env python


import sys
import cPickle as pickle

from ..tools.gtTools import *

from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationSite:
    def __init__(self,genePos,geneName,geneStrand,gtfCnts,TYPE):

        self.genePos    =  genePos
        self.geneName   =  geneName
        self.geneStrand =  geneStrand
        self.type=TYPE
    

        self.refBase = gtfCnts[0]
        self.cov     = float(sum(gtfCnts[1::]))
        myCnts = gtfCnts[1::]
        

        maxDiff = 0
        for x in range(len(myCnts)):
            if x == baseSwap(self.refBase):
                self.refRate = myCnts[x] / self.cov
            else:
                if myCnts[x] > maxDiff:
                    maxDiff=myCnts[x]; self.mutRate = myCnts[x] / self.cov;  self.mutBase = baseSwap(x);   self.mutCnts = myCnts[x]


    def findGenomicLocation(self,spliceTable):
        k=0
        
        while self.genePos > spliceTable[k][0][1]: k+=1

        if self.genePos < spliceTable[k][0][0]:

            if self.type == 'INTRONIC': self.valid=True
                

            minDists=sorted([[abs(spliceTable[k][0][0]-self.genePos),'EX-START',spliceTable[k][0][0],spliceTable[k][1][0]],[abs(self.genePos-spliceTable[k-1][0][1]),'Ex-END',spliceTable[k-1][0][1],spliceTable[k-1][1][1]]])[0]
            
            #print minDists
            #if self.geneStrand     == '-': self.hgPos = spliceTable[k][1][1] + (spliceTable[k][0][1] - self.genePos) #+1
            #else:                          self.hgPos = spliceTable[k][1][1] - (spliceTable[k][0][1] - self.genePos) #+1

        else:
            if self.type != 'EXONIC':
                self.valid =False
                return
            else:
                self.valid = True
            minDists=sorted([[abs(s[0]-self.genePos),'EX-START',s[0],s[1]] for s in spliceTable[k][2]] + [[abs(s[0]-self.genePos),'EX-END',s[0],s[1]] for s in spliceTable[k][3]])[0]


            #if self.geneStrand     == '-': self.hgPos = spliceTable[k][1][1] + (spliceTable[k][0][1] - self.genePos) #+1
            #else:                          self.hgPos = spliceTable[k][1][1] - (spliceTable[k][0][1] - self.genePos) #+1

        if self.geneStrand     == '-': self.hgPos = spliceTable[k][1][1] + (spliceTable[k][0][1] - self.genePos) #+1
        else:                          self.hgPos = spliceTable[k][1][1] - (spliceTable[k][0][1] - self.genePos) #+1
        
        self.spliceDist = minDists[0]; self.spliceType = minDists[1];  self.spliceSite = (minDists[2],minDists[3])
            

