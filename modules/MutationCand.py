#!/usr/bin/env python


import sys
import cPickle as pickle
#from MutationSite import *
from Sequence  import *
from Utilities import *

#from ..tools.gtTools import *

from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationCand:
    def __init__(self,pos):

       
        self.pos = pos
        self.cnts = [0,0,0,0]
        self.diffs = 0
        self.match = 0
        self.refs  = []
        self.quals = [[],[],[],[]]; self.idxs  = [[],[],[],[]]
        self.seqType = None
        self.mapType  = [0,0]
        self.annoType = [0,0]
        self.codonOffset = None

    def record(self,rType,aType,readIdx,readBase,refBase,qualScr):
        if readBase != refBase: self.diffs+=1
        bNum = baseSwap(readBase)
        self.cnts[bNum]+=1; self.quals[bNum].append(ord(qualScr)-37); self.idxs[bNum].append(readIdx)

        if rType == "UNIQ":     self.mapType[0]+=1
        else:                   self.mapType[1]+=1
        if aType == "SINGLE":   self.annoType[0]+=1
        else:                   self.annoType[1]+=1


    def addRef(self,ref):

        if self.seqType == "EXONIC":
            self.refBase = ref[2]
            self.codonSeq = ref[2-self.codonOffset: (2-self.codonOffset)+3]
            self.codonSeq[self.codonOffset] = self.codonSeq[self.codonOffset].lower()
        else:
            self.refBase = ref
            self.codonSeq = "NNN"


    def evalStats(self):

        self.cov = sum(self.cnts); topCnt = 0
        for i in range(len(self.cnts)):

            if baseSwap(i) == self.refBase:
                self.refCnt = self.cnts[i]
                self.refIdx = i
            else:
                if self.cnts[i] > topCnt:
                    topCnt  = self.cnts[i]
                    self.altCnt = self.cnts[i];
                    self.altIdx = i
         
        self.altBase = baseSwap(self.altIdx)
        self.noiseRate   = (self.cov - (self.altCnt + self.refCnt)) / float (self.cov)
        self.altRate = self.altCnt / float(self.cov)
        self.refRate = self.refCnt / float(self.cov)

        if self.refCnt == 0:
            self.urnRef =0; self.qualRef =0
        else:
            self.urnRef = self.urnScr(self.idxs[self.refIdx]); self.qualRef = sum(self.quals[self.refIdx]) / float(self.refCnt)

        self.urnAlt = self.urnScr(self.idxs[self.altIdx]); self.qualAlt = sum(self.quals[self.altIdx]) / float(self.altCnt)



    def evalJxns(self,jData):
        self.spliceClass = "MAJ"
        if self.seqType == "EXONIC":
            
            startDist = self.pos - jData[0][0]; endDist = jData[0][1] - self.pos; self.spliceClass = "MAJ"; self.start_tuple = (jData[0][0],jData[1][0]); self.end_tuple = (jData[0][1],jData[1][1])
            
            if jData[1][0] < jData[1][1]:   self.hgPos = jData[1][0] + startDist
            else:                           self.hgPos = jData[1][0] - startDist

            for s in jData[2]:
                if s[0] <= self.pos and self.pos - s[0] < startDist:     startDist = self.pos - s[0]; self.start_tuple = s;  self.spliceClass = "MIN"
            
            for s in jData[3]:
                if s[0] >= self.pos and s[0] - self.pos < endDist:       endDist = s[0] - self.pos; self.end_tuple = s;      self.spliceClass = "MIN"

            if startDist < endDist:     self.spliceDist = startDist; self.spliceType = "BEG"
            else:                       self.spliceDist = endDist; self.spliceType  = "END"
        
        elif self.seqType == "INTRONIC":

            self.spliceDist = self.pos - jData[0][0]; self.spliceType = "END"
            self.end_tuple = (jData[0][0],jData[1][0])
 
            if jData[1][0] < jData[1][1]:       self.hgPos = jData[1][0] + self.spliceDist
            else:                               self.hgPos = jData[1][0] - self.spliceDist

            if jData[0][1] - self.pos < self.spliceDist:
                self.spliceType = "BEG";    self.spliceDist = jData[0][1] - self.pos;   self.start_tuple = (jData[0][1],jData[1][1])

        elif self.seqType == '3P-FLNK':
            self.spliceDist = self.pos - jData[1] + 1
            if jData[0] == "+": self.hgPos = jData[2] + self.spliceDist
            else:               self.hgPos = jData[2] - self.spliceDist
            self.spliceType = "END"
            self.end_tuple = (jData[1]-1,jData[2])

        elif self.seqType == '5P-FLNK':
            self.spliceDist = int(fabs(self.pos))
            if jData[0] == "+": self.hgPos = jData[2] - self.spliceDist
            else:               self.hgPos = jData[2] + self.spliceDist
            self.spliceType = "BEG"
            self.start_tuple = (0,jData[2])

            #########################################################

            
         









    @staticmethod
    def urnScr(idxs):
        X=len(idxs); Y=len(set(idxs)); Y+=1
        if X<=Y: return 0.999
        if X < 5:
            if X == Y+1: return 0.05
            else:        return 0.01
        if X <= 10:
            if X==Y+1:
                if X==5: return 0.1
                if X==6: return 0.15
                if X==7: return 0.20
                if X==8: return 0.25
                if X==9: return 0.30
                if X==10: return 0.40
            else:
                return 0.045
        else:
            Y+=1
            if X==Y: return 0.999
            elif X==Y+1:
                if X>20: return 0.999
                else: return (Y-10)*0.05 + 0.35
            elif X==Y+2:
                if X>25: return 0.99
                else: return (Y-10) *0.05 + 0.05 
            elif X==Y+3:
                if X > 30: return 0.99
                else:   return (Y-10) * (0.05) + 0.01
            else:
                return 0.5
                
             

           
            





        









