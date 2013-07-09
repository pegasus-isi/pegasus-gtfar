#!/usr/bin/env python


import sys
from SamRead import *
import cPickle as pickle
from collections import defaultdict as dd


##########################################################################################################################################
#####################################################  SAMFILE-FILE CLASS START ##########################################################
##########################################################################################################################################

class SamFile:
    def __init__(self,fileHandle,prefix,strandSpecific=None): #,rlen,lowQ,avgQ,minTrim):
        self.fname = open(fileHandle)
        self.prefix = prefix
        
        
        self.nextRead = self.fname.readline().split()
        
        
        self.open = True
        self.index = 0
        

        if strandSpecific=='16' or strandSpecific=='-':
            self.reqStrand='16'
        elif strandSpecific=='0' or strandSpecific=='+':
            self.reqStrand='0'
        else:
            self.reqStrand=None
        
        while self.nextRead[0][0]=="@":
            self.nextRead=self.fname.readline().split()
        
        self.minOverlap=8
        self.showAMBIG=True
        self.read=''



##############################################################################################################
############################################   ADD KEY   #####################################################
##############################################################################################################

    def addKey(self,keyFile,EXONIC,INTRONIC):
        if fileExtension(keyFile)=="pickle":
            self.key=pickle.load(open(keyFile,'rb'))
        else:
            self.key={}
            if EXONIC:
                for line in open(keyFile):
                    line=line.split()
                    refType=line[0].split("|")[6][0:4]
                    if refType=="FLNK" or refType=="ITRN":
                        continue
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]
            if INTRONIC:
                for line in open(keyFile):
                    line=line.split()
                    refType=line[0].split("|")[6][0:4]
                    if refType=="ITRN" or refType=="GENE":
                        for i in range(1,len(line)):
                            line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                        self.key[line[0]]=[line[1],line[2],line[3]]
                    elif refType=="FLNK":
                        for i in range(1,len(line)):
                            line[i]=[[x.split("-")[0],x.split("-")[1]] for x in line[i].split(",")]
                            for j in range(len(line[i])):
                                for k in range(len(line[i][j])):
                                    if line[i][j][k][0]=="N":
                                        line[i][j][k]=int(line[i][j][k][1::])*-1
                                    else:
                                        line[i][j][k]=int(line[i][j][k])
                        self.key[line[0]]=[line[1],line[2],line[3]]
                    else:
                        continue


##############################################################################################################
############################################   ADD READ   ####################################################
##############################################################################################################
            
    def getRead(self):
        firstID=self.nextRead[0]; multiReads=[]
        while self.nextRead[0] == firstID:
            multiReads.append(self.nextRead)
            self.nextRead = self.fname.readline().split()
            if self.nextRead==[]:
                self.open=False
                break
        ### STRAND SPECIFICITY ###

        SENSE=None
        if self.reqStrand:
            tmpStrands=set([])
            for m in multiReads:
                tmpStrands.add(m[1])
            if len(tmpStrands)>1:
                m=0; SENSE=True
                while m<len(multiReads):
                    if multiReads[m][1]!=self.reqStrand:
                        multiReads.remove(multiReads[m])
                    else:
                        m+=1
            elif tmpStrands.pop() != self.reqStrand:
                SENSE=False
            else:
                SENSE=True
        self.index+=1
        location_keys=[]
        for i in range(len(multiReads)):
            s=multiReads[i][2].split("|")
            sType= "|".join(s[0:len(s)-1])+"|"
            location_keys.append(self.key[sType])
        self.read = SamRead(multiReads,location_keys,SENSE)
        

        
##############################################################################################################
####################################   STORING READ DATA  ####################################################
##############################################################################################################

    
    def storeGeneExpression(self):
        if self.index==1:
            self.cntDict={}; self.geneSeen=set([])
        
        for i in range(len(self.read.locations)):
            for mapping in self.read.locations[i][1::]:
                geneID=mapping[0]
                if geneID not in self.cntDict:
                    self.cntDict[geneID]=[0,0,0];    self.geneSeen.add(geneID)

                if len(self.read.locations)>1:            self.cntDict[geneID][2]+=1
                elif len(self.read.locations[i][1::])>1:  self.cntDict[geneID][1]+=1
                else:                                self.cntDict[geneID][0]+=1

    def storeFeatureExpression(self):
        if self.index==1:
            self.featureCnt={}

        for i in range(len(self.read.data)):
            featureName="|".join(self.read.data[i])
            if featureName not in self.featureCnt:
                self.featureCnt[featureName]=[0,0,0]

            if len(self.read.locations)>1: self.featureCnt[featureName][2]+=1

            elif len(self.read.data)>1: self.featureCnt[featureName][1]+=1

            else:  self.featureCnt[featureName][0]+=1


    def storeSpliceLocations(self):
        if self.index==1:
            self.spliceFeatures=dd(int)
    
        if self.read.junctionSpan:
            self.spliceFeatures[self.read.junctionSpan]+=1


##############################################################################################################
####################################   PROCESS READ  #########################################################
##############################################################################################################


    def processRead(self):
        self.read.relocate()
        self.read.gatherJxns(self.minOverlap)
        self.read.visualizeStrings()
        
        if self.read.sense != False:
            self.storeGeneExpression()
            self.storeFeatureExpression()
            self.storeSpliceLocations()


##############################################  WRITING STUFF ########################################


    def writeToBam(self):
        if self.index==1:
            self.samOut=open(self.prefix+'_viz.sam','w')

        if len(self.read.samStrings)==1:

            self.samOut.write('%s\n' % self.read.samStrings[0])

        elif self.showAMBIG:

            for s in self.read.samStrings:

                self.samOut.write('%s\n' % s)


    def writeGeneLocations(self):
        
        if self.index==1:
            self.readOut=open(self.prefix+'_gene.locs','w')
            if self.reqStrand != None:
                self.antiOut = open(self.prefix+'_antiSense.locs','w')
        
        if self.read.sense != False:

            for s in self.read.locStrings:
                self.readOut.write('%s\n' % s)

        else:
            for s in self.read.locStrings:
                self.antiOut.write('%s\n' % s)


    def printSplices(self,spliceOut):
        for s in self.spliceFeatures.keys():
            spliceOut.write("%s %s\n" % (s,self.spliceFeatures[s]))


    def printFeatures(self,featureOut):
        for s in self.featureCnt.keys():
            featureOut.write("%s %s %s %s\n" % (s,self.featureCnt[s][0],self.featureCnt[s][1],self.featureCnt[s][2]))

    def printRPKMs(self,rpkmOut):

        for key in self.key.keys():
            key=key.split("|")
            if key[6]=="GENE":
                k="|".join(key[0:6])
                exSize=int(key[7]); gSize=int(key[8])
                if k in self.geneSeen:
                    rpkmOut.write("%s %s %s | %s %s %s\n" % (k,exSize,gSize,self.cntDict[k][0],self.cntDict[k][1],self.cntDict[k][2]))
                else:
                    rpkmOut.write("%s %s %s | %s %s %s\n" % (k,exSize,gSize,0,0,0))


