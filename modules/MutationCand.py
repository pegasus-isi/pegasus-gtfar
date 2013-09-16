#!/usr/bin/env python


import sys
from Sequence  import *
from Utilities import *


from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationCand:
    def __init__(self,chr,pos,ref,mutCnts,minCov,minRate):

        self.minCov, self.minRate = minCov,minRate 
        
        self.chr, self.pos, self.ref = chr,pos,ref
        self.alt = baseSwap(mutCnts.index(max(mutCnts)))
        self.refScrs = [0,0,dd(int)]
        self.altScrs = [0,0,dd(int)]
        self.noiseCnt = 0.0
        self.genePositions  = dd(int)
        self.spliceFeatures = dd(lambda: [0,0,0,0])
        


    def update(self,readPos,read,qual,geneLoc,geneFeature):
        idx = self.pos - readPos;
        readBase,qualScr = read[idx],qual[idx]
        readEmbed = min(idx,len(read)-idx)
        if readBase == self.ref:    self.refScrs[0] +=1; self.refScrs[1] +=ord(qualScr); self.refScrs[2][readEmbed]+=1
        elif readBase == self.alt:  self.altScrs[0] +=1; self.altScrs[1] +=ord(qualScr); self.altScrs[2][readEmbed]+=1
        else:                       self.noiseCnt   += 1

       

        mygene = geneLoc.split(":")

        self.genePositions[mygene[0]] = int(mygene[1].split(",")[0]) + idx
        
        self.spliceFeatures[geneFeature][baseSwap(readBase)]+=1
    
    
    def evaluate(self,specHandle,chrNum):
      
        self.valid = False

        self.cov = self.refScrs[0] + self.altScrs[0] + self.noiseCnt

        self.altRate = smartDivide(self.altScrs[0],self.cov); self.refRate = smartDivide(self.refScrs[0],self.cov) ; self.noiseRate = smartDivide(self.noiseCnt,self.cov)

        self.altQual = smartDivide(self.altScrs[1],self.altScrs[0]); self.refQual = smartDivide(self.refScrs[1],self.refScrs[0])
        
        if self.cov < self.minCov or self.altRate < self.minRate or (self.refRate > 0.1 and self.noiseRate > 0.1):   return False

        self.altUrn = self.urnProb(len(self.altScrs[2]),sum(self.altScrs[2].values()))
        self.refUrn = self.urnProb(len(self.refScrs[2]),sum(self.refScrs[2].values()))

        if self.altUrn < 0.1 or (self.altQual < self.refQual/2.0): return False

        self.valid = True

        geneSpots=" ".join([":".join([g,str(self.genePositions[g])]) for g in self.genePositions.keys()])
        spliceFeats = " ".join([f for f in self.spliceFeatures])
        
        print self.chr,chrNum,self.pos,self.ref,self.alt,self.cov,self.altRate,self.altUrn,"GENES:",geneSpots,"FEATS:",spliceFeats
        
        spliceTmp = []

        if len(self.spliceFeatures) > 1:
            altIndex = baseSwap(self.alt)
            for f in self.spliceFeatures:
                spliceCnts = self.spliceFeatures[f]
                spliceCov  = sum(spliceCnts)
                if spliceCov > 5:
                    tmpRate = float(spliceCnts[altIndex]) / spliceCov
                    spliceTmp.append([tmpRate,spliceCnts,f])
            if len(spliceTmp) > 1:
                spliceTmp.sort()
                if spliceTmp[0][0] + 0.2 < spliceTmp[-1][0]:
                    specHandle.write("%s %s %s %s CANDIDATE: " % (self.chr,self.pos,self.ref,self.alt))
                    for s in spliceTmp:
                        specHandle.write("%s %s %s %s " % (s[2],s[0],sum(s[1]),",".join([str(x) for x in s[1]])))
                    specHandle.write("\n")


    @staticmethod
    def urnProb(X,Y):  # = X = number of buckets; Y = number of observations 

        X+=1
        if X > 25: return 0.99
        if X > 20:
            if Y < 50: return 0.99
            if Y < 60: return 0.65
            if Y < 70: return 0.35
            if Y < 80: return 0.10
            return 0.01
        if X > 15:
            if Y < 30: return 0.99
            if Y < 40: return 0.70
            if Y < 50: return 0.55
            if Y < 60: return 0.45
            return 0.01
        if X > 10:
            if Y < 20:  return 0.99
            if Y < 30:  return 0.75
            if Y < 40:  return 0.25
            return 0.01
        if X>5:
            if X+2 >= Y: return 0.99
            if X+3 >= Y: return 0.80
            if X+4 >= Y: return 0.60
            if X+5 >= Y: return 0.20
            return 0.01
        if X+1 >= Y:     return 0.99
        if X+2 >= Y:     return 0.50
        if X+3 >= Y:     return 0.05
        return 0.01

                
             

           
            





        







             

           
            





        









