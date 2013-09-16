#!/usr/bin/env python


import sys
from MapLine import *
from MapRead import *
from Sequence import *
from Utilities import *
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs

##########################################################################################################################################
#####################################################  MAPDATA CLASS START      ##########################################################
########print#############################################################################################################################

class MapData:

    def __init__(self,prefix,keyFile,refType,fName):
 

        self.prefix, self.refType,self.mapFile = prefix,refType,fName
        self.keyHandle = open(keyFile)
        
        self.spliceExp, self.chrCnt, self.typeCnt = dd(int), dd(int), dd(int)
       
        self.invalidCnt, self.readCnt, self.senseCnt,self.spliceCnt,self.uniqCnt = 0,0,0,0,0
        
        self.geneKey, self.biasTable  = dd(lambda: ["MULTI",[0,0,0,0]]), dd(lambda:[0 for i in xrange(100)])

        if self.refType != "INTERGENIC":
            for line in self.keyHandle:
                keyLine    = line.split()[0].split("|")
                if keyLine[6][0:4] == "EXNS":
                    self.geneKey[keyLine[0]] = ["|".join(keyLine[0:5]),[0,0,0,0]]
                
            
            self.samOut    = open(self.prefix+'_vis.sam','w')
            self.geneOut   = open(self.prefix+'_gene.cnts','w') 
            self.spliceOut = open(self.prefix+'_splice.cnts','w')
            self.process     = self.processGeneRead
            self.printResult = self.printGeneResult
        else:
            self.process = self.processGenomeRead
            self.printResult = self.printGenomeResult

##########################################################################################################################################
#####################################################  INIT FUNCTIONS           ##########################################################
##########################################################################################################################################
 

    def processGeneRead(self,mapRead):
        
        if mapRead.invalid:
            self.invalidCnt +=1 
            return 
        self.readCnt +=1; exprIdx = 2
        if mapRead.sense:
            exprIdx =0; self.senseCnt +=1    
        if mapRead.geneUniq:
            self.geneKey[mapRead.geneLoc[0]][1][exprIdx]+=1
        else:
            self.geneKey[",".join(sorted(mapRead.multiGenes))][1][exprIdx+1]+=1
        
        if mapRead.hgUniq:
            self.samOut.write("%s\n" % mapRead.samString())
            
            if mapRead.geneUniq and mapRead.sense:
                mapRead.jxnSearch()

                if mapRead.spliced:
                    self.spliceCnt +=1
                    for s in mapRead.spliceJxns: self.spliceExp[s]+=1
                for read in mapRead.splitReads:
                    sys.stdout.write("%s\n" % " ".join(read))
        if self.refType == "EXONIC" and mapRead.fivePrimeBias:
            self.biasTable[mapRead.geneLoc[0]][int((mapRead.fivePrimeBias[0]*100.0) / mapRead.fivePrimeBias[1])]+=1








    def processGenomeRead(self,mapRead):
        if mapRead.invalid:
            self.invalidCnt +=1 
            return
        self.readCnt +=1
        self.senseCnt += 1

        if mapRead.hgUniq:
            mapRead.hgJxnSearch()
            if mapRead.spliced:
                self.spliceCnt+=1
                for s in mapRead.spliceJxns: self.spliceExp[s]+=1
            for read in mapRead.splitReads:
                sys.stdout.write("%s\n" % " ".join(read))

            #print "YO",mapRead.hgLoc,mapRead.geneLoc,mapRead.read


















    def printGeneResult(self):

    ### --- GENE EXPRESSION ### 

        for gene in self.geneKey:
            geneData,geneCnts = self.geneKey[gene]; geneSum = sum(geneCnts)
            if geneData == "MULTI":
                tmpData = "/".join([self.geneKey[name][0] for name in gene.split(",")])
                self.geneOut.write("%s %s %s\n" % (tmpData,self.refType," ".join([str(s) for s in geneCnts])))
            else:
                self.geneOut.write("%s %s %s\n" % (geneData,self.refType," ".join([str(s) for s in geneCnts])))
                if geneSum > 0:
                    self.chrCnt[geneData.split("|")[3]]  += geneSum
                    self.typeCnt[geneData.split("|")[4]] += geneSum
 
        for s in self.spliceExp:    self.spliceOut.write("%s %s %s\n" % (s,self.refType,self.spliceExp[s]))
        
        if self.refType == "EXONIC" and len(self.biasTable) > 1:
            self.biasOut = open(self.prefix+'_bias.graph','w')
            for k in self.biasTable.keys():
                kCnts = self.biasTable[k]
                if sum(kCnts)>1000:
                    for i in range(1,100):
                        self.biasOut.write("%s %s %s\n" % (k,i,kCnts[i]))


        self.statsOut = open(self.prefix+'.stats','w')
        self.statsOut.write("QC-Statistics: %s \n" % self.mapFile)
        self.statsOut.write("Total-Reads[valid/invalid] %s %s \n" % (self.readCnt, self.invalidCnt))
        self.statsOut.write("Valid-Reads[sense/antisense] %s %s \n" % (self.senseCnt , self.readCnt-self.senseCnt))
        self.statsOut.write("Valid-Reads[spliced/unspliced] %s %s \n" % (self.spliceCnt , self.readCnt-self.spliceCnt))
        self.statsOut.write("MitoChondrial-Reads %s \n" % self.chrCnt["chrM"])
        
        for t in self.typeCnt:  self.statsOut.write("%s-Reads %s\n" % (t.upper(),self.typeCnt[t]))

         
            

    def printGenomeResult(self):


        if self.spliceCnt > 0:
            self.spliceOut = open(self.prefix+'_splice.cnts','w')
            for s in self.spliceExp:
                self.spliceOut.write("%s %s %s\n" % (s,self.refType,self.spliceExp[s]))

        self.statsOut = open(self.prefix+'.stats','w')
        self.statsOut.write("QC-Statistics: %s \n" % self.mapFile)
        self.statsOut.write("Total-Reads[valid/invalid] %s %s \n" % (self.readCnt, self.invalidCnt))
        self.statsOut.write("Valid-Reads[sense/antisense] %s %s \n" % (self.senseCnt , self.readCnt-self.senseCnt))
        self.statsOut.write("Valid-Reads[spliced/unspliced] %s %s \n" % (self.spliceCnt , self.readCnt-self.spliceCnt))
        self.statsOut.write("MitoChondrial-Reads %s \n" % self.chrCnt["chrM"])
        

         
            


##############################################################################################################
####################################   PRINTING DATA #########################################################
##############################################################################################################

    






































































                
            


    def close(self):

        self.mapHandle.close()


        if self.samOut: self.samOut.close()
        if self.mapOut: self.mapOut.close()
        if self.geneMaps > 0: self.geneOut.close()
        if self.spliceMaps > 0: self.spliceOut.close()
        #if self.novelMaps > 0: self.novelOut.close()






