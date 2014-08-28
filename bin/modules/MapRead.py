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

    def __init__(self,fileHandle,fileType,strandSpecific,tag):


        self.mapMax = 50

        if strandSpecific == "OPPOSITE":    self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False,("+","."): True, ("-","."): True, ("X","+"): True, ("X","-"): True}
        elif strandSpecific == "MATCH":     self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True,("+","."): True,("-","."): True, ("X","+"): True, ("X","-"): True}
        else:                               self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True,("+","."): True,("-","."): True, ("X","+"): True, ("X","-"): True}

    
        self.samStrand = {'0':'+','16':'-','+':'0','-':'16'}


        self.handle, self.tag= fileHandle, tag

        self.minOVERLAP=5

        
        self.fileOpen = True
        


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
            elif len(self.rawLine[2].split(":")[0].split("|"))==1:
                if len(self.rawLine[5].split("M"))==2 and len(self.samHeader)>1:
                    for h in self.samHeader: print h
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
            samStartData.extend(['GT:i:'+GENOME_UNIQUE,'TT:i:'+GENE_UNIQUE,CLASS,'GN:i:None','AN:i:None','FM:i:INTERGENIC','SN:i:True','RT:i:'+self.tag])
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
        self.disambiguateMaps()
        self.loadOutputStrings()

    def loadOutputStrings(self):
      
        tmpOutput=[]
        for m in self.mapData:
            self.readPrint, self.qualPrint = self.readSeq,self.qual
            if len(m[1])>2:
                if m[1][1]-m[1][0]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readPrint,self.qualPrint=self.readPrint[(m[1][1]-m[1][0])+1::],self.qualPrint[(m[1][1]-m[1][0])+1::]
                    else:
                        self.readPrint,self.qualPrint=self.readPrint[0:len(self.readPrint)-(m[1][1]-m[1][0]+1)],self.qualPrint[0:len(self.readPrint)-(m[1][1]-m[1][0]+1)]
                    m = (m[0],m[1][2::],m[2],m[3],m[4])
                if m[1][-1]-m[1][-2]<self.minOVERLAP:
                    if m[2]=="+":
                        self.readPrint,self.qualPrint=self.readPrint[0:len(self.readPrint)-(m[1][-1]-m[1][-2]+1)],self.qualPrint[0:len(self.readPrint)-(m[1][-1]-m[1][-2]+1)]
                    else:
                        self.readPrint,self.qualPrint=self.readPrint[(m[1][-1]-m[1][-2])+1::],self.qualPrint[(m[1][-1]-m[1][-2])+1::]
                    m=(m[0],m[1][0:len(m[1])-2],m[2],m[3],m[4])
            
            
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
                    sys.exit()

            elif len(m[1])>2:   m[0][5]=m[0][5]+":JXN="+"-".join([str(m[1][i-1])+','+str(m[1][i]) for i in range(2,len(m[1])-1,2)])
                

            elif len(m[1])==2 and (m[0][5] == "KJXN" or m[0][5] =="NJXN"):
                EDGE=False
                for i in range(len(m[4])):
                    pair = m[4][i]
                    if m[1][0] >= pair[0] and m[1][1] <= pair[1]:
                        EDGE=True
                        m[0][5]="EXON:COORDS="+"-".join([str(s) for s in pair])
                       # else:   m[0][5]="EXON:COORDS="+"-".join([str(s) for s in pair])
                        break
                if not EDGE:
                    print m,self.readID
                    print "WTF EXON HUH"
                    sys.exit()
            elif m[0][5] != "FILTER":
                print "UNKNOWN THING",m[0][5]
                sys.exit()


            cigar="".join([str(m[1][i]-m[1][i-1]+1)+"M" if i%2==1 else str(m[1][i]-m[1][i-1]-1)+"N" for i in xrange(1,len(m[1]))])
            samStrand=self.samStrand[m[2]]
            if m[2]=="-":
                self.readPrint = reverse_complement(self.readPrint)
                self.qualPrint = self.qualPrint[::-1]
            
            if (m[0][4]==m[2]) != m[3]:
                if samStrand == '0': samStrand='16'
                elif samStrand == '16': samStrand='0'
                else:
                    print "WTF"
                    sys.exit()

            
            startOutput = [self.readID,samStrand,m[0][3],m[1][0],'255',cigar,"*",0,0,self.readPrint,self.qualPrint,'NM:i:'+str(self.subs)]
            endOutput   = ['CL:i:'+m[0][5],'GN:i:'+m[0][0],'AN:i:'+m[0][1],'FM:i:'+m[0][2],'SN:i:'+str(m[3]),'RT:i:'+self.tag]
            tmpOutput.append([startOutput,endOutput])
        if len(tmpOutput) == 1:
            uniqOutput = tmpOutput
        else:
            tmpOutput.sort()
            uniqOutput = [tmpOutput[0]]
            for i in range(1,len(tmpOutput)):
                if tmpOutput[i][0]!=uniqOutput[-1][0] or tmpOutput[i][1][1:3] != tmpOutput[-1][1][1:3]:
                    uniqOutput.append(tmpOutput[i])

        if len(uniqOutput) == 1:
            self.outputStrings = [uniqOutput[0][0]+["GT:i:1.0,0","TT:i:1.0,0"]+uniqOutput[0][1]]
        else:
            self.outputStrings = []
            genomicDict = dd(int); featureDict = dd(int)
            for s in uniqOutput:
                genomicDict[(s[0][2],s[0][3])]+=1
                featureDict[(s[1][1],s[1][2])]+=1 
            for s in uniqOutput:
                gCnt = genomicDict[(s[0][2],s[0][3])]; fCnt = featureDict[(s[1][1],s[1][2])]
                if gCnt == 1 and fCnt == 1:
                    self.outputStrings.append(s[0]+["GT:i:0,1.0","TT:i:0,1.0"]+s[1])
                elif gCnt == len(uniqOutput):
                    self.outputStrings.append(s[0]+["GT:i:"+str(1.0/len(uniqOutput))+",0","TT:i:0,"+str(1.0/fCnt)]+s[1])
                elif fCnt == len(uniqOutput):
                    self.outputStrings.append(s[0]+["GT:i:1.0,"+str(1.0/gCnt)   , "TT:i:"+str(1.0/len(uniqOutput))+",0"]+s[1])
                else:
                    self.outputStrings.append(s[0]+["GT:i:1.0,"+str(1.0/gCnt)   , "TT:i:0,"+str(1.0/fCnt)]+s[1])
                    



            


    def printFeatureData(self):
        
        if len(self.outputStrings) < self.mapMax:
            for S in self.outputStrings:
                print "\t".join([str(s) for s in S])

                



    def disambiguateMaps(self):

        self.multiGenes = None

        # MapData Description
        # MapData is a list with 5 elements
        # MapData[0] = Feature information: [ GeneID, AltID, GeneDescription, chromosome, geneStrand, featureType ]
        # MapData[1] = Map Location in genome coordinates 
        # MapData[2] = Map Strand (+ or -) on genome 
        # MapData[3] = Sense Information (True or False) - Always True if Data is not strand specific
        # MapData[4] = Feature Location on Genome ( for example the start and stop position of the exon )


        if len(self.mapData)==1: return 
        #0) REMOVE DUPLICATES
            
        #a) sort data and remove elements which are not identical to previous


        self.mapData.sort()
        self.mapData = self.mapData[0:1]+[self.mapData[i] for i in range(1,len(self.mapData)) if self.mapData[i]!=self.mapData[i-1]];
        if len(self.mapData)==1: return 
        
        
        #1) CHECK FOR MULTI SENSE 
        #a) If more the one value exists in the sense field, keep only those which are in the sense orientation 
        
        if len(set([m[3] for m in self.mapData]))>1:
            self.mapData=[m for m in self.mapData if m[3]]
            if len(self.mapData)==1: return 
        
        
        #2) Prioritize Mapping Features ( Filter is preferential to Exonic Alignment which is preferential to intergenic, etc)
        
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
        
        
        if len(set([tuple(m[0]) for m in self.mapData])) == 1 and len(set([m[1] for m in self.mapData]))>1:
            # If all mapping locations have the same feature data #
            if max([len(m[1]) for m in self.mapData])>2 and len(set([(m[1][1:len(m[1])-1]) for m in self.mapData])) == 1:
                # If all mapping locations contain the same interior locations and more than two locations # 
                tmp_starts=(max([m[1][0] for m in self.mapData]),min([m[1][0] for m in self.mapData])); tmp_ends=(max([m[1][-1] for m in self.mapData]),min([m[1][-1] for m in self.mapData]))
                if tmp_starts[0]-tmp_starts[1] == tmp_ends[0]-tmp_ends[1]:
                    tmp_offset = tmp_starts[0]-tmp_starts[1]
                    if tmp_offset*2 < len(self.readSeq)/2.0:
                        self.mapData = [(self.mapData[0][0],(tmp_starts[0],)+self.mapData[0][1][1:len(self.mapData[0][1])-1]+(tmp_ends[1],),self.mapData[0][2],self.mapData[0][3],self.mapData[0][4])]
                        self.readSeq,self.qual = self.readSeq[tmp_offset:len(self.readSeq)],self.qual[tmp_offset:len(self.readSeq)]
            elif len(set([m[1][-1] for m in self.mapData])) == 1:
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
            elif len(set([m[1][0] for m in self.mapData])) == 1:
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

