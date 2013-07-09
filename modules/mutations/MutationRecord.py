#!/usr/bin/env python


import sys
import cPickle as pickle
from MutationSite import *
from ..tools.gtTools import *

from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationRecord:
    def __init__(self,fileHandle,prefix,TYPE,coverage,diffRate): 
        self.fname = open(fileHandle)
        self.line = self.fname.readline().split()
        self.open = True
        self.chr = self.line[3]; self.gene= self.line[5]; self.index = 0
        
        self.fType = TYPE; self.prefix = prefix
        self.mutFile = None; self.geneCandidates=[]
        self.coverage_parameter = coverage;   self.rate_parameter = diffRate;     self.endLength=-1
    def nextGene(self):
        if self.line == []:
            self.gene = "NULL"; self.chr = "NULL"
        else:
            while self.line[5]== self.gene:
                self.line = self.fname.readline().split()
            self.chr = self.line[3]
            self.gene= self.line[5]
            self.geneCandidates = []; self.flankCandidates=[]
            self.index +=1 

    def findCands(self,gtfRecord):
        self.gtfGene = gtfRecord.geneKey[self.gene];    self.gtfGene.getSeqFromChr(gtfRecord.seq,0);    gtfCnts = [[self.gtfGene.seq[n],0,0,0,0,0] for n in range(len(self.gtfGene.seq))]
        while self.line[5] == self.gene:
            mapLocs=[int(s) for s in self.line[6].split(",")]
            k=0; n=0; x=0; seqLen=len(self.line[8])
           
            if self.fType == 'EXONIC' or (self.fType=='INTRONIC' and mapLocs[0]>=0 and mapLocs[-1]+seqLen< self.gtfGene.length):
                
                for base in self.line[8]:
                    x+=1

                    ### NOTICE WHAT THE CONDITIONS REQUIRE (NOT IN FIRST 3 OR LAST THREE BASES OF READ, ALSO NOT AT THE SPLICE SITE )

                    if x>=self.endLength and x<=(seqLen-self.endLength):   #and n > 0 and mapLocs[k]+n < mapLocs[k+1] +100:
                        gtfCnts[mapLocs[k]+n][baseSwap(base)+1]+=1
                    if mapLocs[k]+n == mapLocs[k+1]:
                        k+=2; n=-1
                    n+=1
            self.line = self.fname.readline().split()
            if self.line == []:
                break 
            
        for i in range(len(gtfCnts)):
            if sum(gtfCnts[i][1::]) > gtfCnts[i][baseSwap(gtfCnts[i][0])+1] + 0:
                self.geneCandidates.append(MutationSite(i,self.gene,self.gtfGene.strand,gtfCnts[i],self.fType))
           
    def evalAndPrint(self):
        snpOut = sys.stdout
        self.gtfGene.findSplicingInfo()

        if self.mutFile == None:
            self.mutFile=open(self.prefix+'.mutations','w')

        for site in self.geneCandidates:
            
            site.findGenomicLocation(self.gtfGene.spliceInfo)

            if site.valid and site.cov > self.coverage_parameter and site.mutRate > self.rate_parameter:
                self.mutFile.write('%s:%s %s %s %s %s %s %s cov: %s rate: %s' % (self.chr,site.hgPos,site.geneStrand,site.type,site.geneName ,site.genePos+1,site.refBase,site.mutBase,site.cov,site.mutRate))
                self.mutFile.write(' | Splice: %s %s %s:%s %s:%s\n' % (site.spliceDist,site.spliceType,site.geneName,site.spliceSite[0],self.chr,site.spliceSite[1]))
                self.gtfGene.seq[site.genePos]=site.mutBase




        


    

