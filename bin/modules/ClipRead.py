#!/usr/bin/env python

import os 
import sys
import re
import difflib
from MapRead import *
#from Sequence import *
#from Utilities import *
from ToolSet import errorQuit
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs


#from ..gtTools.seq import split *


##########################################################################################################################################
#####################################################  MAPPED READ  CLASS START ############################################################
##########################################################################################################################################
##########################################################################################################################################
#####################################################  MAPLINE-FILE CLASS START ##########################################################
########print##################################################################################################################################

class ClipReads:

    def __init__(self,fileHandle,fileType,tag):


        self.mapMax = 50

    

        self.handle, self.tag= fileHandle, tag

        self.minOVERLAP=5

        self.splice_dict= dd(lambda: dd(int))
        self.fileOpen = True
        
        self.spliceDict = dd(int)

        #myLine = self.handle.readline().split()
        self.rawLine = self.handle.readline().split()
        myLine=self.rawLine
        self.samHeader=[]
        if len(self.rawLine) == 0:
            self.fileOpen = False
            self.format = None
            return
        if self.rawLine[0][0] == "@":
            while self.rawLine[0][0] == "@":
                self.samHeader.append("\t".join(self.rawLine))
                #print "\t".join(self.rawLine)
                self.rawLine = self.handle.readline().split()


        self.format = fileType
        if self.format  == "sam":
            self.format = "SAM"
            self.seqLen = len(self.rawLine[9])
            if len(self.rawLine[2].split(":")[0].split("|"))==6:
                self.getNextRead = self.getSamGeneRead
                self.printData   = self.printGeneData
            elif len(self.rawLine[2].split(":")[0].split("|"))==1:
                self.getNextRead = self.getSamGenomeRead
                self.printData   = self.printGenomeData
        else:
            errorQuit(".mapping or .sam extension required")


        
##############################################################################################################
############################################   READ RELOCATION  ##############################################
##############################################################################################################

    

    def relocate(self,locs,fStart):
        read_remain = len(self.readSeq)-1
        offset = 0; outLoc = ()
        pDist = 0 
        for pair in locs:
            pDist = pair[1]-pair[0]+1
            if pDist <= fStart:
                fStart -= pDist
            elif pDist > fStart + read_remain:
                return outLoc + (pair[0]+fStart,pair[0]+fStart+read_remain)
            else:
                outLoc += (pair[0]+fStart,pair[1])
                read_remain -= (pDist-fStart) 
                fStart = 0

    def cigarRelocate(self,locs,fStart,cigar):
        ## TEMPORARY CODE ###
        if len(locs)==1:
            cLocs=[[int(c.split("N")[0]),int(c.split("N")[1])] for c in ("0N"+cigar[0:-1]).split("M")]
            return tuple([item for pair in [[(fStart+locs[0][0])+cLocs[i][0]+(sum([sum(cLocs[j]) for j in range(i)]))-1,(fStart+locs[0][0])+(sum([sum(cLocs[j]) for j in range(i+1)]))-2] for i in range(len(cLocs))] for item in pair])
    
    def genomeRelocate(self,fStart,cigar):
        cLocs=[[int(c.split("N")[0]),int(c.split("N")[1])] for c in ("0N"+cigar[0:-1]).split("M")]
        return tuple([item for pair in [[fStart+cLocs[i][0]+(sum([sum(cLocs[j]) for j in range(i)]))-1,fStart+(sum([sum(cLocs[j]) for j in range(i+1)]))-2] for i in range(len(cLocs))] for item in pair])
         

    

    def getSamGenomeRead(self):
        if len(self.rawLine)==12: self.readID,self.readSeq,self.subs,self.qual,self.mapData=self.rawLine[0],self.rawLine[9],int(self.rawLine[11].split(":")[-1]),self.rawLine[10],[]
        elif len(self.rawLine)==13: self.readID,self.readSeq,self.subs,self.qual,self.mapData=self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[]
        else:
            print "WTF"; sys.exit()
        samData=self.rawLine[0:len(self.rawLine)-1]
        while self.rawLine[0] == self.readID:
            self.mapData.append((self.rawLine[1],self.rawLine[2],self.rawLine[3],self.rawLine[4],self.rawLine[5],self.genomeRelocate(int(self.rawLine[3]),self.rawLine[5]),True,True))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break

        if len(self.mapData)>1:
            self.mapData=[]
            self.outputStrings = []
        else:
            mData = self.mapData[0]
            self.splice_dict[mData[1]][(mData[5][1],mData[5][2])]+=1



    def getSamGeneRead(self):
        self.readID,self.readSeq,self.subs,self.qual,self.mapData,self.mapNum  = self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[],1
        while self.rawLine[0] == self.readID:
            feature=self.rawLine[2].split(":")
            data=self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5])
            
            self.mapData.append((feature[0].split("|"), self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5]),True,True))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break
        
        if len(self.mapData)>1:
            self.mapData = [] 
            self.outputStrings = []
            sys.exit()
        else:
            mData = self.mapData[0]
            spot_data = tuple([s for s in mData[0]])
            self.splice_dict[spot_data][(mData[1][1],mData[1][2])]+=1
        return 



    

    def printGeneData(self):
        for g in self.splice_dict:
            ID1='\"'+g[0]+'\";'
            ID2='\"'+g[1]+'\";'
            TYPE='\"'+g[2]+'\";'
            my_chr = g[3]
            for s in self.splice_dict[g]:
                KIND='\"'+str(self.splice_dict[g][s])+'\";'
                outputString = [my_chr,"GTFAR","novel-jxn-cand",s[0],s[1],".",g[4],".","gene_id",ID1,"transcript_id",ID1,"gene_type",TYPE,"gene_status",KIND,"gene_name",ID2,"transcript_type",TYPE,
                    "transcript_status",KIND,"transcript_name",ID2]

                print "\t".join([str(s) for s in outputString])

    def printGenomeData(self): 
        ID1='\"'+'INTERGENIC'+'\";'
        ID2='\"'+'INTERGENIC'+'\";'
        TYPE='\"''INTERGENIC''\";'

        for c in self.splice_dict:
            for s in self.splice_dict[c]:
                KIND='\"'+str(self.splice_dict[c][s])+'\";'
                outputString = [c,"GTFAR","novel-jxn-cand",s[0],s[1],".",'+',".","gene_id",ID1,"transcript_id",ID1,"gene_type",TYPE,"gene_status",KIND,"gene_name",ID2,"transcript_type",TYPE,
                    "transcript_status",KIND,"transcript_name",ID2]
                print "\t".join([str(s) for s in outputString])
                



