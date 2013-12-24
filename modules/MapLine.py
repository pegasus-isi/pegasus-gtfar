#!/usr/bin/env python


import sys
import re
from MapRead import *
from Sequence import *
from Utilities import *
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs

##########################################################################################################################################
#####################################################  MAPLINE-FILE CLASS START ##########################################################
########print##################################################################################################################################

class MapLines:

    def __init__(self,mapFile):
 
        self.fileName = mapFile
        ##
        try:
            self.handle = open(self.fileName)
        except IOError:
            errorQuit('ERROR: File '+self.fileName+' does not exist')
        
        self.open = True
        
        myLine = self.handle.readline().split()
        if len(myLine) == 0:
            self.open = False
            self.format = None
            return
            #errorQuit("Warning: "+self.fileName+' is empty')
        if myLine[0][0] == "@":
            while myLine[0][0] == "@": myLine = self.handle.readline().split()



        if self.fileName.split(".")[-1] == "mapping":    self.format = "MAPPING"
        elif self.fileName.split(".")[-1] == "sam":      self.format = "SAM"
        else:                                           errorQuit(".mapping or .sam extension required")

        
        
        
        if self.format == "MAPPING":

            self.seqLen = len(myLine[1]) - 1
            self.rawLine = myLine
            self.next = self.nextMappingLine
            self.next()


        if self.format == "SAM":
            self.refData = myLine[2].split("|")
            if len(self.refData) == 1:
                self.rName,self.samStrand,self.fName,self.fPos,self.cigar,self.read,self.qual,self.subs = myLine[0],myLine[1],myLine[2],int(myLine[3]),myLine[5],myLine[9],myLine[10],int(myLine[11].split(":")[-1])
                self.geneID,self.refStrand,self.chr,self.seqLen = None,"+",self.fName,len(self.read)-1
                self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
                self.next = self.nextGenomeLine

            else:
                self.rName,self.samStrand,self.fName,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = myLine
                self.fPos = int(self.fPos)-1; self.subs= int(self.subs.split(":")[-1]); self.refData = self.fName.split("|")
                self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
                self.seqLen = len(self.read) - 1 
                self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
                self.next = self.nextSamLine

        
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



    def nextMappingLine(self):
        
        try:
            self.rName,self.seq,refFeature,fPos,ref,strand,self.subs,self.locs,self.qual = self.rawLine
            refFeature = refFeature.split(":")
            #geneID,geneAlt,geneDesc,chr,geneStrand,refDesc = refFeature[0].split("|")
            geneID,geneAlt,geneGroup,refChr,refStrand,refType = refFeature[0].split("|")
            
            #genegeneID,geneAlt,
            #geneLoc = (refDesc,geneAlt,geneID,geneStrand,geneDesc)
            geneLoc = (geneAlt,geneID,geneGroup,refChr,refStrand,refType)
            
            refLoc  = (refChr,strand,self.relocate([ [int(r.split('-')[0]),int(r.split('-')[1])] for r in refFeature[1].split("|")],int(fPos)))
            #geneLoc = (refDesc,geneAlt,geneID,geneStrand,geneDesc)
            self.data = (refLoc,geneLoc,ref)




            self.rawLine = self.handle.readline().split()

        except ValueError:
            self.open = False; self.rName = None
         
    def nextSamLine(self):
            
        try:
            self.rName,self.samStrand,self.fName,self.fPos,self.blank,self.cigar,self.blank,self.blank,self.blank,self.read,self.qual,self.subs,self.blank = self.handle.readline().split()
            self.fPos = int(self.fPos) - 1; self.subs= int(self.subs.split(":")[-1]); self.refData = self.fName.split("|")
            self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
            self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
        except ValueError:
            self.open = False; self.rName = None

    def nextGenomeLine(self):
        myLine = self.handle.readline().split()
        try:
            self.rName,self.samStrand,self.fName,self.fPos,self.cigar,self.read,self.qual,self.subs = myLine[0],myLine[1],myLine[2],int(myLine[3]),myLine[5],myLine[9],myLine[10],int(myLine[11].split(":")[-1])
            self.geneID,self.refStrand,self.chr,self.seqLen = None,"+",self.fName,len(self.read)-1
            self.fLocs = cigarToLoc(int(self.fPos),self.cigar)
        except IndexError:
            self.open = False; self.rName = None


'''
class InfoTemplate:
    def __init__(self,TYPE,data):

        if TYPE == "REFLOC":
            self.chr,self.strand,self.locs = data

        if TYPE == "GENEMAP":
            self.type,self.id,self.hugo,self.strand  = data




'''











 






























