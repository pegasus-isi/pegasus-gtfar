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

class MapReads:

    def __init__(self,fileHandle,fileType,strandSpecific):
 
        if strandSpecific == "OPPOSITE":    self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False}
        elif strandSpecific == "MATCH":         self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True}
        else:                               self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True}

    
        self.samStrand = {'0':'+','16':'-','+':'0','-':'16'}


        self.handle = fileHandle; self.gapped = False

        self.minOVERLAP=5


        #try:                self.handle = open(self.fileName)
        #except IOError:     errorQuit('ERROR: File '+self.fileName+' does not exist')
        
        self.open = True

        self.fileOpen = True
        


        #myLine = self.handle.readline().split()
        self.rawLine = self.handle.readline().split()
        myLine=self.rawLine
        self.samHeader=[]
        if len(self.rawLine) == 0:
            self.open = False
            self.fileOpen = False
            self.format = None
            return
        if self.rawLine[0][0] == "@":
            while self.rawLine[0][0] == "@":
                self.samHeader.append("\t".join(self.rawLine))
                #print "\t".join(self.rawLine)
                self.rawLine = self.handle.readline().split()


        self.format = fileType
        if self.format  == "mapping":
            self.seqLen = len(self.rawLine[1]) - 1
            self.format = "MAPPING"
            if len(self.rawLine[2].split(":")[0].split("|"))==6:
                self.getNextRead = self.getMappingFeatureRead
                self.printData   = self.printFeatureData
            else:
                errorQuit("No Mapped Reads")
        
        elif self.format  == "sam":
            self.format = "SAM"
            self.seqLen = len(self.rawLine[9])
            if len(self.rawLine[2].split(":")[0].split("|"))==6:
                self.getNextRead = self.getSamFeatureRead
                self.printData   = self.printFeatureData
            elif self.rawLine[2][0:3]=="chr":
                if len(self.rawLine[5].split("M"))==2 and len(self.samHeader)>1:
                    for h in self.samHeader: print h
                self.getNextRead = self.getSamGenomeRead
                self.printData   = self.printGenomeData
        else:
            errorQuit(".mapping or .sam extension required")


            ### FIX THIS AND MAKE IT MORE CONCISE ###!!!!
            #self.readID,self.readSeq,self.subs,self.qual  = self.rawLine[0],self.rawLine[1],self.rawLine[6],self.rawLine[8]
            
#        if self.format == "SAM":
#            print "YO"
#            cigarSpot = myLine[5].split("N"); self.refData = myLine[2].split("|")
#            if len(cigarSpot) > 1 and len(self.refData)>1:
#                self.style = 'gapFeature'
#                self.seqlen = len(myLine[9]); self.rawLine = myLine; self.next = self.nextGappedLine; self.next()
                #print myLine 

#            elif len(self.refData) == 1:
#                self.rName,self.samStrand,self.fName,self.fPos,self.cigar,self.read,self.qual,self.subs = myLine[0],myLine[1],myLine[2],int(myLine[3]),myLine[5],myLine[9],myLine[10],int(myLine[11].split(":")[-1])
#                self.geneID,self.refStrand,self.chr,self.seqLen = None,"+",self.fName,len(self.read)-1
#                self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
#                self.next = self.nextGenomeLine

     
 #           else:
 #               self.rName,self.samStrand,self.fName,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = myLine
 #               self.fPos = int(self.fPos)-1; self.subs= int(self.subs.split(":")[-1]); self.refData = self.fName.split("|")
 #               self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
 #               self.seqLen = len(self.read) - 1 
 #               self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
 #               self.next = self.nextSamLine

        
##############################################################################################################
############################################   READ RELOCATION  ##############################################
##############################################################################################################

    

    def relocate(self,locs,fStart):
        read_remain = self.seqLen
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
        #print "\t".join(self.rawLine)
        if len(self.rawLine)==12: self.readID,self.readSeq,self.subs,self.qual,self.genomeData=self.rawLine[0],self.rawLine[9],int(self.rawLine[11].split(":")[-1]),self.rawLine[10],[]
        elif len(self.rawLine)==13: self.readID,self.readSeq,self.subs,self.qual,self.genomeData=self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[]
        else:
            print "WTF"; sys.exit()
        samData=self.rawLine[0:len(self.rawLine)-1]
        while self.rawLine[0] == self.readID:
            self.genomeData.append((self.rawLine[1],self.rawLine[2],self.rawLine[3],self.rawLine[4],self.rawLine[5],self.genomeRelocate(int(self.rawLine[3]),self.rawLine[5])))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break
    def printGenomeData(self):

        GENE_UNIQUE="0.0,0.0"

        if len(self.genomeData)==1:  GENOME_UNIQUE="1.0,0.0"
        else:                        GENOME_UNIQUE="0.0,"+str(1.0/len(self.genomeData))

        
        for g in self.genomeData:
            if len(g[5])==2:
                CLASS="CL:i:INTERGENIC:REG="+str(g[5][0])+","+str(g[5][1])
            else:
                CLASS="CL:i:INTERGENIC:JXN="+str(g[5][1])+","+str(g[5][2])
            samStartData= [self.readID,g[0],g[1],g[2],g[3],g[4],"*",0,0,self.readSeq,self.qual,'NM:i:'+str(self.subs)]
            samStartData.extend(['GT:i:'+GENOME_UNIQUE,'TT:i:'+GENE_UNIQUE,CLASS,'GN:i:None','AN:i:None','FM:i:INTERGENIC','SN:i:True'])
            print "\t".join([str(s) for s in samStartData])


    def getSamFeatureRead(self):
        self.readID,self.readSeq,self.subs,self.qual,self.mapData,self.mapNum  = self.rawLine[0],self.rawLine[9],0,self.rawLine[10],[],1
        while self.rawLine[0] == self.readID:
            feature=self.rawLine[2].split(":")
            data=self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5])
            
            self.mapData.append((feature[0].split("|"), self.cigarRelocate([[int(x) for x in r.split("-")] for r in  feature[1].split("|")],int(self.rawLine[3]),self.rawLine[5]),self.samStrand[self.rawLine[1]],True))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break

        if len(self.mapData)>1:
            print "HI"
            sys.exit()
        return 


    def getMappingFeatureRead(self):
        self.readID,self.readSeq,self.subs,self.qual,self.mapData,self.mapNum  = self.rawLine[0],self.rawLine[1],self.rawLine[6],self.rawLine[8],[],1

        while self.rawLine[0] == self.readID:
            feature=self.rawLine[2].split(":")
            tmpLocs=[[int(x) for x in r.split("-")] for r in feature[1].split("|")]
            self.mapData.append((feature[0].split("|"), self.relocate( tmpLocs, int(self.rawLine[3])),self.rawLine[5],self.strandKey[self.rawLine[5],feature[0].split("|")[4]],tmpLocs))
            self.rawLine = self.handle.readline().split()
            if len(self.rawLine)==0:
                self.fileOpen=False
                break
        #self.mapType = "UNIQUE"
        if len(self.mapData)>1: self.disambiguateMaps()
       


    def printFeatureData(self):
        
        ### RATE GENE/GENOME UNIQUENESS ###

        if len(set([tuple(m[0][0:3]) for m in self.mapData]))==1: GENE_UNIQUE=str(1.0/len(self.mapData))+",0.0"
        else:                                                     GENE_UNIQUE="0.0"+str(1.0/len(self.mapData))
        
        if len(set([m[1] for m in self.mapData]))==1:             GENOME_UNIQUE=str(1.0)+",0.0"
        else:                                                     GENOME_UNIQUE="0.0,"+str(1.0/len(self.mapData))

        for m in self.mapData:

            ### FIX SMALL SPLICING OFFSETS ###
            if len(m[1])>2:
                if m[1][1]-m[1][0]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readSeq,self.qual=self.readSeq[(m[1][1]-m[1][0])+1::],self.qual[(m[1][1]-m[1][0])+1::]
                    else:
                        self.readSeq,self.qual=self.readSeq[0:len(self.readSeq)-(m[1][1]-m[1][0]+1)],self.qual[0:len(self.readSeq)-(m[1][1]-m[1][0]+1)]
                    m = (m[0],m[1][2::],m[2],m[3],m[4])
                if m[1][-1]-m[1][-2]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readSeq,self.qual=self.readSeq[0:len(self.readSeq)-(m[1][-1]-m[1][-2]+1)],self.qual[0:len(self.readSeq)-(m[1][-1]-m[1][-2]+1)]
                    else:
                        self.readSeq,self.qual=self.readSeq[(m[1][-1]-m[1][-2])+1::],self.qual[(m[1][-1]-m[1][-2])+1::]
                    m=(m[0],m[1][0:len(m[1])-2],m[2],m[3],m[4])
            

            SPLICESTR='None'

            
            if m[0][5]=="EXON":
                m[0][5]="EXON:COORDS="+"-".join([str(s) for s in m[4][0]])
            elif m[0][5]=="ITRN" or m[0][5] == "INTRON":
                if m[1][0] >= m[4][1][0] and m[1][1] <= m[4][1][1]: m[0][5]="ITRN:COORDS="+"-".join([str(s) for s in m[4][1]])
                elif m[1][0] >= m[4][1][0]: m[0][5]="IJXN:EX-INT="+"-".join([str(s) for s in m[1]])
                elif m[1][1] <= m[4][0][1]: m[0][5]="IJXN:INT-EX="+"-".join([str(s) for s in m[1]])
                elif m[1][0] <= m[4][1][0] and m[1][1] >= m[4][1][0]:   m[0][5]="IJXN:EX-INT="+"-".join([str(m[4][0][1]),str(m[4][1][0])])
                elif m[1][0] <= m[4][1][1] and m[1][1] >= m[4][2][0]:   m[0][5]="IJXN:INT-EX="+"-".join([str(m[4][1][1]),str(m[4][2][0])])
                else:
                    print "WTF"
                    print m
                    print ""
                    sys.exit()

            elif len(m[1])>2:   m[0][5]=m[0][5]+":JXN="+"-".join([str(m[1][i-1])+','+str(m[1][i]) for i in range(2,len(m[1])-1,2)])
                

            elif len(m[1])==2 and (m[0][5] == "KJXN" or m[0][5] =="NJXN"):
                EDGE=False
                for i in range(len(m[4])):
                    pair = m[4][i]
                    if m[1][0] >= pair[0] and m[1][1] <= pair[1]:
                        EDGE=True
                        if i == 0:  m[0][5]="EXON:END="+str(pair[1])
                        elif i == len(m[4])-1: m[0][5]="EXON:START="+str(pair[0])
                        else:   m[0][5]="EXON:COORDS="+"-".join([str(s) for s in pair])
                        break
                if not EDGE:
                    print m,self.readID
                    print "WTF EXON HUH"
                    sys.exit()
            elif m[0][5] != "FILTER":
                print m
                print "UNKNOWN THING",m[0][5]
                sys.exit()


            
            cigar="".join([str(m[1][i]-m[1][i-1]+1)+"M" if i%2==1 else str(m[1][i]-m[1][i-1]-1)+"N" for i in xrange(1,len(m[1]))])
            
            samStartData= [self.readID,self.samStrand[m[2]],m[0][3],m[1][0],'255',cigar,"*",0,0,self.readSeq,self.qual,'NM:i:'+str(self.subs)]

 

            #samStartData.extend(['GT:i:'+GENOME_UNIQUE,'SP:i:'+SPLICESTR,'TT:i:'+GENE_UNIQUE,'CL:i:'+m[0][5],'GN:i:'+m[0][0],'AN:i:'+m[0][1],'FM:i:'+m[0][2],'SN:i:'+str(m[3])])
            samStartData.extend(['GT:i:'+GENOME_UNIQUE,'TT:i:'+GENE_UNIQUE,'CL:i:'+m[0][5],'GN:i:'+m[0][0],'AN:i:'+m[0][1],'FM:i:'+m[0][2],'SN:i:'+str(m[3])])
            
            print "\t".join([str(s) for s in samStartData])



                



    def disambiguateMaps(self):

        ### 0) REMOVE DUPLICATES
        self.mapData.sort();    self.mapData = self.mapData[0:1]+[self.mapData[i] for i in range(1,len(self.mapData)) if self.mapData[i]!=self.mapData[i-1]];
        if len(self.mapData)==1: return
        ### 1) CHECK FOR MULTI SENSE ### 
        if len(set([m[3] for m in self.mapData]))>1:   self.mapData=[m for m in self.mapData if m[3]]
        ### 2) Prioritize Map Codes ###
        mapCodes = [m[0][-1] for m in self.mapData]
        if "FILTER" in mapCodes:
            self.mapData=[m for m in self.mapData if m[0][-1]=="FILTER"]
            if len(self.mapData)>1:
                tmpData=[m for m in self.mapData if m[0][1] != "CHR"]
                if len(tmpData)>1:  tmpData=[t for t in tmpData if t[0][1][0:6] != "MT-RNR"]
                if len(tmpData)>0: self.mapData=tmpData

        elif "EXON" in mapCodes: self.mapData=[m for m in self.mapData if m[0][-1]=="EXON"]
        elif "KJXN" in mapCodes: self.mapData=[m for m in self.mapData if m[0][-1]=="KJXN"]
        if len(self.mapData)==1: return 
        ### 3) Merge Multi-Genes ###
        if len(set([m[1] for m in self.mapData]))==1 and len(set([m[0][3] for m in self.mapData]))==1:
            self.mapData=[([",".join(list(set([self.mapData[i][0][j] for i in range(len(self.mapData))]))) for j in range(len(self.mapData[0][0]))], self.mapData[0][1],self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
            if len(self.mapData)==1: return

        ### 4) Merge Mutlti-Spots --- ###
        
        if len(set([tuple(m[0]) for m in self.mapData])) == 1:
            if len(set([(m[1][1:len(m[1])-1]) for m in self.mapData])) == 1:
                tmp_starts=(max([m[1][0] for m in self.mapData]),min([m[1][0] for m in self.mapData])); tmp_ends=(max([m[1][-1] for m in self.mapData]),min([m[1][-1] for m in self.mapData]))
                if tmp_starts[0]-tmp_starts[1] == tmp_ends[0]-tmp_ends[1]:
                    tmp_offset = tmp_starts[0]-tmp_starts[1]
                    if tmp_offset*2 < len(self.readSeq)/2.0:
                        self.mapData = [(self.mapData[0][0],(tmp_starts[0],)+self.mapData[0][1][1:len(self.mapData[0][1])-1]+(tmp_ends[1],),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        self.qual = self.readSeq[tmp_offset:len(self.readSeq)-tmp_offset],self.qual[tmp_offset:len(self.readSeq)-tmp_offset]
            else:
                if len(set([m[1][-1] for m in self.mapData])) == 1:
                    tmpLocs=[]; k=-1;
                    while len(set([m[1][k] for m in self.mapData]))==1:
                        tmpLocs.append(self.mapData[0][1][k])
                        k-=1
                    if len(tmpLocs)%2==1:   tmpLocs.append(max([m[1][k] for m in self.mapData]))
                    
                    tmpLocs.reverse(); trimDist=len(self.readSeq)-sum([(tmpLocs[i]-tmpLocs[i-1])+1 for i in range(1,len(tmpLocs),2)])
                    if trimDist < len(self.readSeq)/2.0:
                        self.mapData=[(self.mapData[0][0],tuple(tmpLocs),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        if self.mapData[0][2] == "-": self.readSeq=self.readSeq[0:len(self.readSeq)-trimDist]; self.qual=self.qual[0:len(self.qual)-trimDist]
                        else:                         self.readSeq=self.readSeq[trimDist::]; self.qual=self.qual[trimDist::]
                        return 
                if len(set([m[1][0] for m in self.mapData])) == 1:
                    tmpLocs = []; k=0; 
                    while len(set([m[1][k] for m in self.mapData]))==1:
                        tmpLocs.append(self.mapData[0][1][k])
                        k+=1

                    if len(tmpLocs)%2==1:   tmpLocs.append(min([m[1][k] for m in self.mapData]))
                    
                    trimDist=len(self.readSeq)-sum([(tmpLocs[i]-tmpLocs[i-1])+1 for i in range(1,len(tmpLocs),2)])


                    if trimDist < len(self.readSeq)/2.0:
                        self.mapData=[(self.mapData[0][0],tuple(tmpLocs),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        
                        
                        
                        if self.mapData[0][2] == "+": self.readSeq=self.readSeq[0:len(self.readSeq)-trimDist]; self.qual=self.qual[0:len(self.qual)-trimDist]
                        else:                         self.readSeq=self.readSeq[trimDist::]; self.qual=self.qual[trimDist::]
                        return
            
                        
                

