#!/usr/bin/env python


import sys
from MutationCand import *
from Sequence import *
from Utilities import *


from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationRecord:
    def __init__(self,mapFile,prefix,coverage,diffRate,SPECIES="HUMAN"):
        self.mapFile = mapFile; self.mapHandle = open(self.mapFile);    self.open = True; self.offset = [0,0]; self.chr = None
        
        self.spliceBuffer = 3; self.coverage_parameter = coverage;   self.rate_parameter = diffRate; self.cands = {}
      

        HUMAN_CHRS = {'chr1': 1, 'chr2': 2, 'chr3': 3,'chr4': 4, 'chr5': 5, 'chr6': 6,'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr10': 10, 'chr11': 11, 'chr12': 12,
                                'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18,'chr19': 19, 'chr20': 20, 'chr21': 21,'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrM': 25}
        
	MONKEY_CHRS = {'chr01': 1, 'chr02a': 2, 'chr02b': 3, 'chr03': 4, 'chr04': 5, 'chr05': 6,'chr06': 7, 'chr07': 8, 'chr08': 9,'chr09': 10, 'chr10': 11, 'chr11': 12, 'chr12': 13,
                    'chr13': 14, 'chr14': 15, 'chr15': 16, 'chr16': 17, 'chr17': 18, 'chr18':19 ,'chr19': 20, 'chrX': 21}


	if SPECIES == "HUMAN": self.chr2Num = dd(lambda: None, HUMAN_CHRS)
	elif SPECIES == "RHESUS" or SPECIES == "MONKEY": self.chr2Num == dd(lambda: None, MONKEY_CHRS)
        
	self.nextLine()
        self.chr = self.mapChr
        if not self.open:       errorQuit("EMPTY LOCATION FILE")
        
	self.specificSnps = open(prefix+'_specific.snps','w')

        #if SPECIES == "HUMAN":
        #    self.chr2Num ={'chr1': 1, 'chr2': 2, 'chr3': 3,'chr4': 4, 'chr5': 5, 'chr6': 6,'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr10': 10, 'chr11': 11, 'chr12': 12,
                              #  'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18,'chr19': 19, 'chr20': 20, 'chr21': 21,'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrM': 25}
        #elif SPECIES == "RHESUS" or SPECIES == "MONKEY":
        #    self.chr2Num = {'chr01': 1, 'chr02a': 2, 'chr02b': 3, 'chr03': 4, 'chr04': 5, 'chr05': 6,'chr06': 7, 'chr07': 8, 'chr08': 9,'chr09': 10, 'chr10': 11, 'chr11': 12, 'chr12': 13,
        #            'chr13': 14, 'chr14': 15, 'chr15': 16, 'chr16': 17, 'chr17': 18, 'chr18':19 ,'chr19': 20, 'chrX': 21}


##########################################################################################################################################################
    
    
    def nextLine(self,SEEK=None):
        if SEEK != None:
            self.mapHandle.seek(SEEK)
            
        tmpOffset = self.mapHandle.tell()
        self.line = self.mapHandle.readline().split()



	while len(self.line) > 0 and self.chr2Num[self.line[1]] == None:
		self.line = self.mapHandle.readline().split()
   
        if len(self.line) > 0:

            self.readID,self.mapChr,self.mapStrand,self.hgStart,self.hgEnd,self.mapRead,self.mapRef,self.mapQual,self.mapSubs,self.geneLoc,self.geneFeature = self.line
            self.hgStart, self.hgEnd, self.mapSubs = int(self.hgStart), int(self.hgEnd), int(self.mapSubs)
                
            if self.mapChr != self.chr:
                if self.chr > self.mapChr:
		    errorQuit("UNSORTED LOCATION FILE")
                self.offset = [self.offset[1],tmpOffset]
	else:
	    self.offset = [ self.offset[1],tmpOffset] 
            self.open = False 
	    self.mapGene,self.mapChr = None, None




###########################################################################################################################################################
           


    def geneCandSearch(self):
        while self.chr == self.mapChr:
            if self.mapSubs > 0 and (len(self.mapRead) > 2*self.spliceBuffer):
                diffs = [(self.hgStart+i,self.mapRef[i],self.mapRead[i]) for i in range(self.spliceBuffer,len(self.mapRead)-self.spliceBuffer) if self.mapRead[i] != self.mapRef[i]]
                for d in diffs:
                    if d[0] not in self.cands.keys():
                        self.cands[d[0]] = [d[1],[0,0,0,0]]
                    self.cands[d[0]][1][baseSwap(d[2])]+=1
            self.nextLine()


    def geneCandCall(self):
        
        candList = []; n = 0
        for c in self.cands:
            if max(self.cands[c][1]) > self.coverage_parameter/2 :    candList.append((c,MutationCand(self.chr,c,self.cands[c][0],self.cands[c][1],self.coverage_parameter,self.rate_parameter  )))
        
        if len(candList) == 0:
            self.chr = self.mapChr
            if self.open: self.mapHandle.seek(self.offset[1])
            return 
        self.cands = {}; candList.sort()
        ##   ------------------------------------------------------------------------  ##
        
	self.nextLine(SEEK = self.offset[0])
	 
        while self.chr == self.mapChr:
            n = 0
            while n < len(candList):
                position,candidate = candList[n]
                if position < self.hgStart:
                    candidate.evaluate(self.specificSnps,self.chr2Num[self.chr])
                    candList.remove(candList[n])
                elif position < self.hgEnd:
                    if self.hgStart + self.spliceBuffer <= position and position <= self.hgEnd - self.spliceBuffer:
                        candidate.update(self.hgStart,self.mapRead,self.mapRef,self.geneLoc,self.geneFeature)
                    n+=1
                else:
                    break
            self.nextLine()
            if len(candList) == 0: break 

        
        for c in candList: c[1].evaluate(self.specificSnps,self.chr2Num[self.chr])
        
        if self.open:
            self.nextLine(SEEK = self.offset[1])
	    self.chr = self.mapChr
        else:
            return
    

