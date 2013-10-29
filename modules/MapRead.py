#!/usr/bin/env python

import os
import sys
import difflib
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs

from Sequence import *

#from ..gtTools.seq import split *


##########################################################################################################################################
#####################################################  MAPPED READ  CLASS START ############################################################
##########################################################################################################################################



class MapRead:
    def __init__(self,mapLines,strandRule):
 
        
        self.strandRule = strandRule
        if strandRule == "OPPOSITE":
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False}
        elif strandRule == "MATCH":
            self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True}
        else:
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True}

       
        self.validRef = dd(lambda: False); self.validRef["chr"] = True
        self.line = mapLines
        self.minOverlap = 5
        

        self.switchStrand = {'0':'+','16':'-','+':'0','-':'16'}
        
        if mapLines.format == "MAPPING":
            self.loadRead = self.pullMappingLine

        
        elif mapLines.format == "SAM":
            if mapLines.refType == "INTERGENIC":
                self.loadRead = self.pullGenomeLine
            else:
                self.loadRead = self.pullSamLine
        

                


    def pullMappingLine(self):

        self.name = self.line.rName
        self.subs = self.line.subs
        self.seq  = self.line.seq
        self.qual = self.line.qual
        
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        self.sense = False 
        self.readLines = [self.line.data]
        #self.readData  = [(self.line.hg,self.line.gene,self.line.ref)]
        self.line.next()
        while self.name == self.line.rName:
            if self.line.data != self.readLines[-1]:
                self.readLines.append(self.line.data)
         #       self.readData.append((self.line.hg,self.line.gene,self.line.ref))
            self.line.next()
        self.processLines()

      
       

    def pullSamLine(self):

        self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.samStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        self.ref = self.read
        readMaps = []
        while self.name == self.line.rName:

            self.relocateTuple(self.key[self.line.fName], self.line.fLocs) 
            readMaps.append([(self.line.chr,self.line.refStrand,self.hgPos),(self.line.geneID,self.genePos),(self.switchStrand[self.line.samStrand],self.line.refStrand)])
            self.line.next()
        if len(readMaps) == 1:
            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]] 
            self.hgUniq, self.geneUniq = True, True 
        else:
            if len(readMaps)<12:
                self.multiMaps = readMaps
                self.disambiguate()
            else:
                self.invalid = True
                return 
        
    
    def pullGenomeLine(self):

        self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.samStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None
        self.ref = self.read
        
        readMaps = []

        while self.name == self.line.rName:
            readMaps.append([(self.line.chr,self.line.refStrand,self.line.fLocs),(self.line.geneID,self.line.fLocs),(self.switchStrand[self.line.samStrand],self.line.refStrand)])
            self.line.next()
        if len(readMaps) == 1:
            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]] 
            self.hgUniq, self.geneUniq = True, True 
        else:
            if len(readMaps)<12:
                self.multiMaps = readMaps
                self.disambiguate()
            else:
                self.invalid = True
                return 
        
    






########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################




    def strandCheck(self):
        if self.strandRule:
            strandScrs = []
            strandScrs = [ self.strandKey[r[0][1],r[1][4]] for r in self.readLines]
            if len(set(strandScrs)) > 1:
                self.readLines = [self.readLines[i] for i in range(len(strandScrs)) if strandScrs[i]]
                return True
            else:
                return strandScrs[0]
        return True

        
    #def removeFalse(self,PRIMARY="EXON",SECONDARY="KJXN"):
    def removeFalse(self,TYPE_FILTERS = ["EXON","KJXN","SPLICE=AG-AT","SPLICE=CT-AC"],GROUP_FILTERS=["exCANON_CANON_exCANON","exCANON_NOVEL_exCANON"]):

        ## FILTER FOR MITO

        
     


        SPECIAL_MAPS = [r for r in self.readLines if r[0][0]=="chrM"]
        if len(SPECIAL_MAPS)>0:
             self.readLines = SPECIAL_MAPS[0:1]
             return 

        ## 0) REMOVE DUPLICATES ##


        self.readLines.sort()
        dupList = [0] + [i for i in range(1,len(self.readLines)) if self.readLines[i] != self.readLines[i-1]]
        if len(dupList) < len(self.readLines):
            self.readLines = [self.readLines[i] for i in dupList]
        
        ## 1) FILTER FOR INVALID CHRS ##
     




        chrVal  = [i for i in range(len(self.readLines)) if self.validRef[self.readLines[i][0][0][0:3]]]
        if len(chrVal) < len(self.readLines) and len(chrVal) > 0:
            self.readLines = [self.readLines[i] for i in chrVal]
            if len(self.readLines) == 1:    return 

        ## 2) FILTER FOR NOVEL/KNOWN JXNS --- VALTYPES HEIRARCHY  ##
        
        types   = [r[1][5] for r in self.readLines]
        if len(set(types))>1:
            for T in TYPE_FILTERS:
                if T in types:
                    valTypes = [i for i in range(len(types)) if types[i] == T]
                    self.readLines = [ self.readLines[i] for i in valTypes]
                    break
                    
        groups = [r[1][2] for r in self.readLines]
        
        if len(set(groups)) >1:
            for G in GROUP_FILTERS:
                if G in groups:
                    valGroups = [i for i in range(len(groups)) if groups[i] == G]
                    self.readLines = [ self.readLines[i] for i in valGroups]
                    break
            
            
        if len(self.readLines) == 1:    return 
       
        ## 2b) FILTER FOR SAME PATH ## 
        if len(set([(r[0][2][0],r[0][2][-1]) for r in self.readLines]))==1 and (self.readLines[0][0][2][0],self.readLines[0][0][2][-1]) in [r[0][2] for r in self.readLines]:
            self.readLines = [[(r[0][0],r[0][1],(r[0][2][0],r[0][2][-1])),r[1],r[2]]  for r in self.readLines]


        ## 3) FILTER FOR POOR NAMING ###

        


        if len(set([r[0] for r in self.readLines])) == 1 and (len(set([r[1][0] for r in self.readLines])) ==1 or len(set([r[1][1] for r in self.readLines]))==1):
                self.readLines = self.readLines[0:1]
                return

        ## 4) FILTER FOR OVERLAPS 
        if len(set([r[1][1] for r in self.readLines])) == 1:
            leftCnt = len(set([r[0][2][0:len(r[0][2])-2] for r in self.readLines])); leftLen= max([r[0][2][-1]-r[0][2][-2] for r in self.readLines])
            rightCnt = len(set([r[0][2][2::] for r in self.readLines])); rightLen= max([r[0][2][1]-r[0][2][0] for r in self.readLines])
            if (leftCnt ==1 and leftLen < self.minOverlap) or (rightCnt == 1 and rightLen < self.minOverlap): #len(set([r[0][2][0:len(r[0][2])-2] for r in DATA])) == 1 and max([r[0][2][-1]-r[0][2][-2]]) <5:
                self.readLines = self.readLines[0:1]
                return



    def trimOffsets(self):
        if not self.hgUniq: return
        if self.locs[1]-self.locs[0]  < self.minOverlap:
            OFFSET=self.locs[1]-self.locs[0]+1
            self.locs = self.locs[2::]
            self.seq=self.seq[OFFSET::]; self.ref=self.ref[OFFSET::]; self.qual = self.qual[OFFSET::]
        if self.locs[-1]-self.locs[-2] < self.minOverlap:
            OFFSET = self.locs[-1]-self.locs[-2]+1
            self.locs = self.locs[0:len(self.locs)-2]
            self.seq = self.seq[0:len(self.seq)-OFFSET]; self.ref = self.ref[0:len(self.ref)-OFFSET]; self.qual = self.qual[0:len(self.qual)-OFFSET]



    @staticmethod 
    def determine_splice_string(geneObj):
        if geneObj.refType == "KJXN":
            SPLICE_TYPE = "exCANON_CANON_exCANON"
        elif geneObj.refType == "NJXN":
            SPLICE_TYPE = "exCANON_NOVEL_exCANON"
        elif geneObj.refType == "ITRN":
            SPLICE_TYPE = "INTRON_EXON_BOUNDARY"
        else:
            SPLICE_TYPE = geneObj.geneGroup
       
        if "exNOVEL" in SPLICE_TYPE.split("_") or "intNOVEL" in SPLICE_TYPE.split("_"):
            SPLICE_DATA = geneObj.refType.split("=")[1]+"|"+SPLICE_TYPE
        else:
            SPLICE_DATA = "CANON"+"|"+SPLICE_TYPE

        return SPLICE_DATA

    def determineSpliceSites(self):
        if not self.hgUniq or len(self.locs)<3:
            self.spliced = False 
        else:
            self.spliced = True
            SPLICE_SITE=",".join([str(self.locs[i-1])+"-"+str(self.locs[i]) for i in range(2,len(self.locs),2)])

            if self.geneUniq:
                self.spliceSite = "|".join([self.chr,self.geneID,self.geneAltID,self.determine_splice_string(self),SPLICE_SITE])
            else:
                myIDS=",".join([m.geneID for m in self.multiGenes])
                myALTS=",".join([m.geneAltID for m in self.multiGenes])
                self.spliceSite = "|".join([self.chr,myIDS,myALTS,self.determine_splice_string(self),SPLICE_SITE])
                    
                




    def processLines(self,SPECIAL="chrM"):
        
        self.hgUniq, self.geneUniq, self.mito, self.ambiguous, self.repetitive,self.spliced = False,False,False,False,False,False
        self.multiLocs, self.multiGenes = None, None
        self.sense = self.strandCheck()

        if len(self.readLines) > 1: self.removeFalse()
       
        self.chr,self.strand,self.locs = self.readLines[0][0]
        self.geneID,self.geneAltID,self.geneGroup,self.refChr,self.refStrand,self.refType = self.readLines[0][1] 
        self.ref = self.readLines[0][2]

        
        if len(set([r[0] for r in self.readLines])) == 1: self.hgUniq   = True
        if len(set([r[1] for r in self.readLines])) == 1: self.geneUniq = True
        

        if not self.hgUniq:
            self.multiLocs =  [InfoTemplate("HGLOC",(r[0],r[1],r[2])) for r in list(set([R[0] for R in self.readLines]))]
            if len(set([r[2] for r in self.readLines])) == 1:  self.repetitive = True
            else:                                              self.ambiguous = True
            
        if not self.geneUniq:
            self.multiGenes = [ InfoTemplate("GENE",(r[0],r[1],r[2],r[3],r[4],r[5])) for r in list(set([R[1] for R in self.readLines]))]
            #print len(self.multiGenes)
        if self.chr == SPECIAL: self.mito = True
        self.trimOffsets()
        self.determineSpliceSites()




        
            


       

        
        


      


#####################################################################################################################################################
    






    def samString(self):
        myChr,myStrand,myPos = self.hgLoc
        cigar = "".join([str(myPos[i]-myPos[i-1]+1)+"M" if i%2==1 else str(myPos[i]-myPos[i-1]-1)+"N" for i in xrange(1,len(myPos))])
        return "\t".join([self.name,self.switchStrand[myStrand],myChr,str(myPos[0]),'255',cigar,'*','0','0',self.read,self.qual,'NM:i:'+str(self.subs)])



















class InfoTemplate:
    def __init__(self,TYPE,data):
        if TYPE == "GENE":
            self.geneID,self.geneAltID,self.geneGroup,self.refChr,self.refStrand,self.refType = data
        
        elif TYPE=="HGLOC":
            self.chr,self.strand,self.locs = data 








































































