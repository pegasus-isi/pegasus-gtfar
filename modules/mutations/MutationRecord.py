#!/usr/bin/env python


import sys
sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')
import cPickle as pickle
from tools.gtTools import *
from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationRecord:
    def __init__(self,fileHandle): #,rlen,lowQ,avgQ,minTrim):
        self.fname = open(fileHandle)
        self.line = self.fname.readline().split()
        self.open = True
        self.chr = self.line[3]
        self.gene= self.line[5]
        self.index = 0
        self.cands=[]
        self.gtfGene=''

    def nextGene(self):
        while self.line[5]== self.gene:
            self.line = self.fname.readline().split()
        self.chr = self.line[3]
        self.gene= self.line[5]
        self.index +=1 

    def call(self,gtfRecord):
        
        
        self.gtfGene = gtfRecord.geneKeys[self.gene]
        self.gtfGene.getSeqFromChr(gtfRecord.seq,1)
        gtfCnts = [[self.gtfGene.seq[n],0,0,0,0,0] for n in range(len(self.gtfGene.seq))]
        while self.line[5] == self.gene:
            mapLocs=[int(s) for s in self.line[6].split("-")]
            k=0; n=0
            for base in self.line[8]:
                if n>=0 and mapLocs[k]+n-1 <= mapLocs[k+1]:
                    gtfCnts[mapLocs[k]+n][baseSwap(base)+1]+=1
                n+=1
                if mapLocs[k]+n-1 == mapLocs[k+1]:
                    k+=2; n=0
            self.line = self.fname.readline().split()
        
        
        for i in range(len(gtfCnts)):
            if sum(gtfCnts[i][1::]) != gtfCnts[i][baseSwap(gtfCnts[i][0])+1]:
                self.cands.append([i,gtfCnts[i]])
           
    def rateAndWrite(self):
        spliceSites=set([])
        for t in self.gtfGene.transcripts:
            for i in range(len(t[1][1])):
                spliceSites.add((tuple(t[1][1][i]),tuple(t[1][2][i])))
        spliceSites=sorted(list(spliceSites))
        print spliceSites[0:5]
        for c in self.cands:
            pos = c[0]
            refBase = c[1][0]
            refSwap = baseSwap(refBase)
            cnts = c[1][1::]
            cov  = float(sum(c[1][1::]))
            mutRate = 0
            for x in range(len(cnts)):
                if x==refSwap:
                    refRate = cnts[x]/cov
                else:
                    if cnts[x]/cov > mutRate:
                        mutRate = cnts[x]/cov
                        mutBase = baseSwap(x)
            print c,"cov",cov
            print "ref",refBase,refRate
            print "mut",mutBase,mutRate



















        
    def addKey(self,keyFile):
        if fileExtension(keyFile)=="pickle":
            self.key=pickle.load(open(keyFile,'rb'))
        else:
            self.key={}
            for line in open(keyFile):
                line=line.split()
                for i in range(1,len(line)):
                    line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                self.key[line[0]]=[line[1],line[2],line[3]]

            
################ NOTE HERE IS THE IDEA FOR GENE EXPRESSION -- ADD THE MULTI GUYS AND THEN GO PRINT THE ZEROS FROM THE KEY, AVOIDING (OR NOTING) ANY THAT HAVE BEEN ADDED AS MULTIS ####


    

