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
        self.handle = open(self.fileName)
        self.open = True
        
        myLine = self.handle.readline().split()
        if len(myLine) == 0:    errorQuit(self.mapFile+' is empty')
        if myLine[0][0] == "@":
            while myLine[0][0] == "@": myLine = self.handle.readline().split()



        if self.fileName.split(".")[-1] == "mapping":    self.format = "MAPPING"
        elif self.fileName.split(".")[-1] == "sam":      self.format = "SAM"
        else:                                           errorQuit(".mapping or .sam extension required")

        refFeature = myLine[2].split("|")

        if len(refFeature) != 1 and len(refFeature) != 7:  errorQuit("BAD REF FEATURE")

        if len(refFeature) == 7:
            fString = refFeature[-1][0:4]
            if   fString  == "ITRN": self.refType = "INTRONIC"
            elif fString  == "EXNS": self.refType = "EXONIC"
            elif fString  == "KJXN": self.refType = "KNOWN-JXNS"
            elif fString  == "NJXN": self.refType = "NOVEL-JXNS"
            elif fString  == "GENE": self.refType = "FULL-GENE"
            else:
                print "WEIRD REF TYPE!"
                errorQuit("BAD REF TYPE")
        else:
            self.refType = "INTERGENIC"


        if self.format == "MAPPING":
            self.rName,self.read,self.fName,self.fPos,self.ref,self.mapStrand,self.subs,self.locs,self.qual = myLine
            self.fPos = int(self.fPos); self.subs = int(self.subs); self.refData = self.fName.split("|")
            self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]

            self.seqLen = len(self.read)-1
            self.next = self.nextMappingLine

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

    
    def nextMappingLine(self):
        
        try:
            self.rName,self.read,self.fName,self.fPos,self.ref,self.mapStrand,self.subs,self.locs,self.qual = self.handle.readline().split()
            self.fPos = int(self.fPos); self.subs = int(self.subs); self.refData = self.fName.split("|")
            self.geneID,self.refStrand,self.chr = self.refData[0],self.refData[2],self.refData[3]
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












































































#################################################################################################################################################################################
########################################################################    PARSING READLINES  ##################################################################################
#################################################################################################################################################################################


    def parseHgLines(self):
        readID = self.mapLines[0][0]
        if self.format == "SAM":
            seq= self.mapLines[0][9]; qual = self.mapLines[0][10]; subs=int(self.mapLines[0][11].split(":")[2])
            mapData = []
            for m in self.mapLines:
                mapLocs=cigarToLoc(int(m[3]),m[5])
                hgTmp = [m[2],mapLocs]
                geneTmp = ["INTERGENIC",mapLocs]
                mapData.append([hgTmp,geneTmp,(self.mappingStrand[m[1]],"+"),[seq,seq],"INTERGENIC"])
            self.read = MapRead(readID,qual,subs,mapData,None,self.refType)
            

         

    def parseSamLines(self):
        readID = self.mapLines[0][0]; seq =  self.mapLines[0][9]; qual = self.mapLines[0][9];  mapData = []; subs=int(self.mapLines[0][11].split(":")[2])
        for i in range(len(self.mapLines)):
            self.featureData = self.mapLines[i][2].split("|")
            chrTmp = self.featureData[3]; nameTmp = self.featureData[0]; featureStrand = self.featureData[2];
            tmpLocs = cigarToLoc(int(self.mapLines[0][3]),self.mapLines[0][5]); geneTmp,hgTmp = [],[]
            myKey = self.key[self.mapLines[0][2]]
            for j in range(0,len(tmpLocs),2):
                geneTuple,hgTuple =self.relocater(myKey,tmpLocs[j],tmpLocs[j+1]-tmpLocs[j])
                geneTmp.extend(geneTuple); hgTmp.extend(hgTuple)
            mapData.append([[chrTmp,hgTmp],[nameTmp,geneTmp],(self.mappingStrand[self.mapLines[i][1]],featureStrand),[seq,seq]])
        self.read = MapRead(readID,qual,subs,mapData,self.strandSpecific,self.refType)


########################################################################################################################################################################################
###########################################################################  END PARSING READLINES  ####################################################################################
########################################################################################################################################################################################


#######################################################################################################################
#######################################################################################################################






    def writeLocations(self):
        samScr = '1'
        self.read.splitLocations()

        if self.refType == "INTERGENIC": return 
        if not self.samOut:
            self.samOut = open(self.prefix+'_vis.sam','w')
        if self.read.valid:
            if self.read.hgUniq and self.read.sense and self.read.geneUniq:
                samScr = '255'
                hg = self.read.hgLocs[0]; samLoc = locToCigar(hg[2]);
                self.samOut.write("%s\n" % "\t".join([self.read.ID,self.samStrand[hg[1][1]],hg[0],str(hg[2][0]),samScr,samLoc,'*','0','0',self.read.seq,self.read.qual,'NM:i:'+str(self.read.subs)]))
            else:
                for i in range(len(self.read.hgLocs)):
                    hg = self.read.hgLocs[i]; samLoc = locToCigar(hg[2])
                    self.samOut.write("%s\n" % "\t".join([self.read.ID,self.samStrand[hg[1][1]],hg[0],str(hg[2][0]),samScr,samLoc,'*','0','0',self.read.seq,self.read.qual,'NM:i:'+str(self.read.subs)]))
        else:
            for i in range(len(self.read.multiGenome)):

                hg =  self.read.multiGenome[i]
                tmpSeq = hg[3][0]; tmpQual = self.read.qual
                samLoc = locToCigar(hg[2])
                if hg[1][0] == "-": tmpQual = tmpQual[-1::-1]
                if len(tmpQual) > len(tmpSeq): tmpQual=tmpQual[0:len(tmpSeq)]
                self.samOut.write("%s\n" % "\t".join([self.read.ID,self.samStrand[hg[1][1]],hg[0],str(hg[2][0]),samScr,samLoc,'*','0','0',tmpSeq,tmpQual,'NM:i:'+str(self.read.subs)]))
                

########################################################################################################################################################################################
#################################################################################  STORING COUNTS   ####################################################################################
########################################################################################################################################################################################

    def storeCnts(self):

        self.read.findSplices()
        for jxn in self.read.spliceTags:
            self.spliceExp[jxn] +=1
        if self.read.sense:     idx = 1
        else:                   idx = 2
        
        self.read.findGeneMaps()
        if len(self.read.geneMaps) > 0:
            self.GENE_EXPRESSION=True
            if len(self.read.geneMaps) == 1:
                self.geneKey[self.read.geneMaps[0]][idx] +=1
            else:
                if not self.read.hgUniq: idx+=1
                self.multiGene[self.read.geneMaps][idx]+=1
           # else:
           #     if not self.read.hgUniq: idx+=2
           #     tmpMap = ",".join([m for m in self.read.geneMaps])
           #     if tmpMap not in self.geneKey:
           #         tmpData = ",".join([self.geneKey[gene][0] for gene in self.read.geneMaps])
           #         self.geneKey[tmpMap] = [tmpData,0,0,0,0]
           #     self.geneKey[tmpMap][idx]+=1

        if self.refType == "INTERGENIC":
            if self.read.valid:
                myChr,roundMaps = self.read.roundChrMaps()
                if myChr not in self.novelAreas:
                    self.novelAreas[myChr] = dd(int)
                for r in roundMaps: self.novelAreas[myChr][r] +=1
        

    def storeQC(self):
        
        self.stats["TOTAL"]+=1

        self.mapSubs[self.read.subs] +=1
        
        if self.read.hgUniq:    self.stats["UNIQUE"]+=1
        else:                   self.stats["AMBIGUOUS"]+=1
        if self.read.spliced:   self.stats["SPLICED"]+=1
        else:                   self.stats["UNSPLICED"]+=1

        if self.strandSpecific != None:
            if self.read.sense: self.stats["SENSE"]+=1
            else:               self.stats["ANTISENSE"]+=1
        if self.refType != "INTERGENIC":
            if self.read.geneUniq:  self.stats["SINGLE-GENE"]+=1
            else:                   self.stats["MULTI-ANNO"]+=1
## GOTTA ADD THIS BACK ###
#       if len(self.read.geneMaps) == 1:
 #           print self.read.geneMaps 
  #          geneData=self.geneKey[self.read.geneMaps[0]][0].split("|")
  #          if geneData[3] == "chrM":
  #              self.mapTypes["MITOCHONDRIAL"]+=1
  #          else:
  #              self.mapTypes[geneData[4]] +=1
        
  #      if self.refType == "EXONICy" and not self.read.spliced and len(self.mapLines)==1 and self.format=="MAPPING":
  #          tmpFeature = self.mapLines[0][2].split("|"); tmpPos = int(self.mapLines[0][3]); fLen = self.key[self.mapLines[0][2]][0][-1][-1]; rLen=len(self.mapLines[0][1])
  #          dist = int(round((float(tmpPos)/(fLen-(rLen-1))) * 100))
  #          if tmpFeature[2] == "-": dist = 100 - dist
  #          fName = "|".join(tmpFeature[0:4])
  #          if fName not in self.bias.keys():
  #              self.bias[fName] = [0 for i in range(101)]
  #          self.bias[fName][dist]+=1







########################################################################################################################################################################################
############################################################################## END STORING COUNTS   ####################################################################################
########################################################################################################################################################################################

########################################################################################################################################################################################
#################################################################################  WRITING COUNTS   ####################################################################################
########################################################################################################################################################################################
    
    
    
    def writeOutCnts(self):
        if len(self.spliceExp)>1:
            self.spliceOut = open(self.prefix+'_splice.cnts','w')
            for s in self.spliceExp.keys():
                self.spliceOut.write("%s %s %s\n" % (s,self.refType,self.spliceExp[s]))
        if self.GENE_EXPRESSION:
            self.geneOut = open(self.prefix+'_gene.cnts','w'); tmpKey = {}
            for gene in self.geneKey:
                geneData = self.geneKey[gene]
                self.geneOut.write("%s | %s uniq: %s %s multi: %s %s [sense/anti]\n" % (geneData[0],self.refType,geneData[1],geneData[2],geneData[3],geneData[4]))
        
        return
        if len(self.novelAreas) > 0:
            self.minNovelCnt = 0; NOVELPASS = False; self.passAreas = {}
            for chr in self.novelAreas:
                self.passAreas[chr]=dd(list)
                for area in self.novelAreas[chr]:
                    if self.novelAreas[chr][area] > self.minNovelCnt:
                        NOVELPASS = True;   self.passAreas[chr][area] = []
            self.novelAreas.clear(); self.mapHandle.close(); self.mapHandle = open(self.mapFile)
            if NOVELPASS:
                self.novelOut = open(self.prefix+'_novelMaps.locs','w')
                for line in self.mapHandle:
                    line=line.split()
                    myChr= line[2]; myPos = int(line[3]); myRound=myPos - myPos%1000; 
                    if self.format == "SAM":    myStrand = self.mappingStrand[line[1]];
                    else:                       myStrand = line[5]
                    if myChr in self.passAreas and myRound in self.passAreas[myChr]:
                        self.passAreas[myChr][myRound].append((myPos,myStrand))
            
                for chr in self.passAreas:
                    for area in self.passAreas[chr]:
                        locs = sorted([x[0] for x in self.passAreas[chr][area]])
                        strands = [x[1] for x in self.passAreas[chr][area]]
                        pos=0;neg=0
                        for s in strands: 
                            if s == "+":
                                pos+=1
                            else:
                                neg+=1
#                        self.novelOut.write('%s %s %s +/- %s %s SPOTS: %s\n' % (chr,area,len(locs),pos,neg," ".join([str(l) for l in locs])))


    def writeQC(self):
        if len(self.bias) > 0:
            self.biasOut = open(self.prefix+'_bias.graph','w')
            for k in self.bias.keys():
                kCnts = self.bias[k]
                if sum(kCnts)>1000:
                    for i in range(1,100):
                        self.biasOut.write("%s %s %s\n" % (k,i,kCnts[i]))
        if len(self.stats) > 0:
            self.statsOut = open(self.prefix+'.stats','w')
            #self.statsOut = sys.stdout
            self.statsOut.write("QC-Statistics for Alignment File ( %s )\n\n" % self.refType)
            self.statsOut.write("Total Reads Processed:     %s\n" % self.stats["TOTAL"])
            subList = ", ".join(sorted([str(k)+" "+str(self.mapSubs[k]) for k in self.mapSubs]))
            self.statsOut.write("Number of Mismatches:      %s\n" % subList)
            
            self.statsOut.write("UnSpliced/Spliced:         %s %s\n" % (self.stats["UNSPLICED"],self.stats["SPLICED"]))
            self.statsOut.write("Genomic Unique/Ambiguous:  %s %s\n" % (self.stats["UNIQUE"],self.stats["AMBIGUOUS"]))
            self.statsOut.write("One/Many Gene Annotations: %s %s\n" % (self.stats["SINGLE-GENE"],self.stats["MULTI-ANNO"]))
            self.statsOut.write("Sense/Antisense:           %s %s\n\n" % (self.stats["SENSE"],self.stats["ANTISENSE"]))
            
        if self.GENE_EXPRESSION:
            typeCnts = sorted([(self.mapTypes[k],k.upper()) for k in self.mapTypes.keys()],reverse=True)
            self.statsOut.write("Alignment Gene Types (ordered by prevalance):\n\n")
            for t in typeCnts:
                self.statsOut.write("%s: %s\n" % (t[1],t[0]))


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






