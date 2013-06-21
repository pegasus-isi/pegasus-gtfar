#!/usr/bin/env python


import sys
from SamRead import *
import cPickle as pickle

##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class SamFile:
    def __init__(self,fileHandle): #,rlen,lowQ,avgQ,minTrim):
        self.fname = open(fileHandle)
        self.nextRead = self.fname.readline().split()
        self.open = True
        self.index = 0
        while self.nextRead[0][0]=="@":
            self.nextRead=self.fname.readline().split()
        
        

        ########## PAD THEM WITH X OR Y BUT NOT SAME DUH #############

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

            
    def storeGeneExp(self,read):
        if self.index==1:
            self.cntDict={}
            self.geneSeen=set([])
        if read.hgUniq and not read.multiAnnos:
            if read.geneLoc[1] not in self.cntDict:
                self.cntDict[read.geneLoc[1]]=[1,0,0]
            else:
                self.cntDict[read.geneLoc[1]][0]+=1
            self.geneSeen.add(read.geneLoc[1])
        elif read.hgUniq and read.multiAnnos:
            tmpAnnos=sorted([anno[1] for anno in read.annos])
            for t in tmpAnnos:
                self.geneSeen.add(t)
                if t not in self.cntDict:
                    self.cntDict[t]=[0,1,0]
                else:
                    self.cntDict[t][1]+=1
        else:
            for k in read.hgKey.keys():
                tmpAnnos=sorted([anno[1] for anno in read.hgKey[k]])
                for t in tmpAnnos:
                    self.geneSeen.add(t)
                    if t not in self.cntDict:
                        self.cntDict[t]=[0,0,1]
                    else:
                        self.cntDict[t][2]+=1

    def printRPKMs(self,rpkmOut):
        for k in self.cntDict.keys():
            rpkmOut.write("%s %s %s %s\n" % (k,self.cntDict[k][0],self.cntDict[k][1],self.cntDict[k][2]))
        for k in self.key.keys():
            tmpKey="|".join(k.split("|")[0:8])
            if tmpKey not in self.geneSeen:
                rpkmOut.write("%s %s %s %s\n" % (tmpKey,0,0,0))





    def getRead(self,strandSpecific=None):
        firstID=self.nextRead[0]; multiReads=[]
        while self.nextRead[0] == firstID:
            multiReads.append(self.nextRead)
            self.nextRead = self.fname.readline().split()
            if self.nextRead==[]:
                self.open=False
                break
        ### STRAND SPECIFICITY ###

        if strandSpecific:
            tmpStrands=set([])
            for m in multiReads:
                tmpStrands.add(m[1])
            if len(tmpStrands)>1:
                m=0
                while m<len(multiReads):
                    if multiReads[m][1]!=strandSpecific:
                        multiReads.remove(multiReads[m])
                    else:
                        m+=1
        self.index+=1
        location_keys=[]
        for i in range(len(multiReads)): location_keys.append(self.key[multiReads[i][2]])
        return SamRead(multiReads,location_keys,self.index)
        
        
################ NOTE HERE IS THE IDEA FOR GENE EXPRESSION -- ADD THE MULTI GUYS AND THEN GO PRINT THE ZEROS FROM THE KEY, AVOIDING (OR NOTING) ANY THAT HAVE BEEN ADDED AS MULTIS ####


    

