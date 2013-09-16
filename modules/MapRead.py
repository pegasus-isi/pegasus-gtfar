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
    def __init__(self,mapLines,readKey,strandRule):
 
        
        self.strandRule = strandRule
        if strandRule == "OPPOSITE":
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): False,("-","-"): False}
        elif strandRule == "MATCH":
            self.strandKey = {("+","-"): False,("-","+"): False,("+","+"): True,("-","-"): True}
        else:
            self.strandKey = {("+","-"): True,("-","+"): True,("+","+"): True,("-","-"): True}

        refAbbrev = {'EXONIC': 'EXNS', 'KNOWN-JXNS': 'KJXN', 'NOVEL-JXNS': 'NJXN', 'INTRONIC': 'ITRN','INTERGENIC': 'GENE','FULL-GENE': 'GENE'}
        
        
        self.line = mapLines
        myType = refAbbrev[mapLines.refType]
        

        self.key = {}
        keyHandle = open(readKey)

        for keyLine in keyHandle:
            keyLine = keyLine.split()
            if keyLine[0].split("|")[6][0:4] == myType:
                self.key[keyLine[0]] = [[[int(x.split('-')[0]),int(x.split('-')[1])] for x in keyLine[i].split(",")] for i in range(1,len(keyLine))]

        self.switchStrand = {'0':'+','16':'-','+':'0','-':'16'}
        
        if mapLines.format == "MAPPING":
            self.loadRead = self.pullMappingLine

        
        elif mapLines.format == "SAM":
            if mapLines.refType == "INTERGENIC":
                self.loadRead = self.pullGenomeLine
            else:
                self.loadRead = self.pullSamLine
        



    def pullMappingLine(self):

        self.name, self.read, self.subs, self.qual, self.firstStrand = self.line.rName, self.line.read, self.line.subs, self.line.qual, self.line.mapStrand
        self.invalid, self.hgUniq, self.geneUniq, self.hgLoc, self.geneLoc, self.fivePrimeBias = False, False, False, None, None, None


        refSet = set([]); readMaps = []
        while self.name == self.line.rName:

            if self.line.mapStrand == "-":     refSet.add(reverse_complement(self.line.ref))
            else:                              refSet.add(self.line.ref)
           
            self.relocate(self.key[self.line.fName],self.line.fPos,self.line.seqLen)
            readMaps.append([(self.line.chr,self.line.refStrand,self.hgPos),(self.line.geneID,self.genePos),(self.line.mapStrand,self.line.refStrand)])
            self.line.next()
        if len(readMaps) == 1:
            self.ref = refSet.pop()
            if self.firstStrand  == "-":
                self.read = reverse_complement(self.read); self.qual = self.qual[-1::]
            
            self.hgLoc, self.geneLoc, self.sense = readMaps[0][0], readMaps[0][1], self.strandKey[readMaps[0][2]]
            self.hgUniq, self.geneUniq = True, True 
        else:
            self.fivePrimeBias = None
            if len(refSet) == 1 and len(readMaps)<12:
                self.ref = refSet.pop()
                if self.firstStrand  == "-":
                    self.read = reverse_complement(self.read); self.qual = self.qual[-1::]
                self.multiMaps = readMaps
                self.disambiguate()
            else:
                self.invalid = True
                self.sense   = True
                return 
        
      
       

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
        
    














      


#####################################################################################################################################################
    






















      


#####################################################################################################################################################
    


    def relocateTuple(self,key,fTuple):

        fKey,gKey,hKey =  key[0],key[1],key[2];  POST=False; gPos = []; hPos =[]

        gGrow = gPos.extend
        hGrow = hPos.extend

        for j in range(0,len(fTuple),2):
            fPos = fTuple[j]; DIST = fTuple[j+1]-fTuple[j]; POST = False
            for i in range(len(fKey)):
                if POST or (fPos >= fKey[i][0] and fPos <= fKey[i][1]):
                    if not POST:
                        myOffset = fPos - fKey[i][0]
                    gGrow([gKey[i][0]+ myOffset, gKey[i][0] + myOffset + DIST  ])

                    if hKey[i][0] < hKey[i][1]:
                        hGrow( [ hKey[i][0] + myOffset, hKey[i][0] + myOffset + DIST ] )
                    else:
                        hGrow( [ hKey[i][0] - myOffset, hKey[i][0] - myOffset - DIST ] )

                    if fKey[i][1] >= fPos + DIST and fKey[i][1] >= fKey[i][0] + DIST:
                        break
                    else:
                        gPos[-1] = gKey[i][1]; hPos[-1] = hKey[i][1]
                        DIST -= (1 + fKey[i][1] - (fKey[i][0] + myOffset) )
                        POST = True; myOffset =0
        self.genePos,self.hgPos = gPos,hPos



########################################################################################################################################################



    def relocate(self,key,fPos,DIST):


        fKey,gKey,hKey =  key[0],key[1],key[2];  POST=False; gPos = []; hPos =[]
        self.fivePrimeBias = (fPos,fKey[-1][-1])

        gGrow = gPos.extend
        hGrow = hPos.extend

        for i in range(len(fKey)):
            if POST or (fPos >= fKey[i][0] and fPos <= fKey[i][1]):
                if not POST:
                    myOffset = fPos - fKey[i][0]
                gGrow([gKey[i][0]+ myOffset, gKey[i][0] + myOffset + DIST  ])

                if hKey[i][0] < hKey[i][1]:
                   hGrow( [ hKey[i][0] + myOffset, hKey[i][0] + myOffset + DIST ] )
                else:
                   hGrow( [ hKey[i][0] - myOffset, hKey[i][0] - myOffset - DIST ] )

                if fKey[i][1] >= fPos + DIST and fKey[i][1] >= fKey[i][0] + DIST:
                    break
                else:
                    gPos[-1] = gKey[i][1]; hPos[-1] = hKey[i][1]
                    DIST -= (1 + fKey[i][1] - (fKey[i][0] + myOffset) )
                    POST = True; myOffset =0
        self.genePos,self.hgPos = gPos,hPos
        #return gPos,hPos



########################################################################################################################################################




    def resolveAntiSense(self):
        tmpMaps = [s for s in self.multiMaps if self.strandKey[s[2]]]

        if len(tmpMaps) > 0:
            self.sense = True
            if len(tmpMaps) == 1:
                self.hgUniq, self.geneUniq = True, True
                self.hgLoc,self.geneLoc= tmpMaps[0][0],tmpMaps[0][1]
            else:
                self.multiMaps = tmpMaps
        else:
            self.sense = False

    def resolveMultipleAnnotations(self):
        self.multiMaps.sort()
        if self.multiMaps[0][0] == self.multiMaps[-1][0]:
            self.hgUniq = True 
            self.hgLoc   = self.multiMaps[0][0]
            self.sense = self.strandKey[self.multiMaps[0][2]]
            self.multiGenes = [myMap[1][0] for myMap in self.multiMaps]
        elif self.multiMaps[0][0][0] == self.multiMaps[-1][0][0]:
            tmpGenes = [myMap[1][0] for myMap in self.multiMaps]
            if len(set(tmpGenes)) == 1:
                self.geneUniq = True
                self.geneLoc = self.multiMaps[0][1]
            else:
                self.multiGenes = tmpGenes
        else:
            self.multiGenes = [myMap[1][0] for myMap in self.multiMaps]
            #self.geneLoc = [geneMap[1] for geneMap in self.multiMaps]
            self.sense = self.strandKey[self.multiMaps[0][2]]
        return


    def disambiguate(self):

        if self.strandRule != None: self.resolveAntiSense()
        else:                       self.sense = True
        if not self.hgUniq:         self.resolveMultipleAnnotations() 
        
        






    def samString(self):

        myChr,myStrand,myPos = self.hgLoc
        cigar = "".join([str(myPos[i]-myPos[i-1]+1)+"M" if i%2==1 else str(myPos[i]-myPos[i-1]-1)+"N" for i in xrange(1,len(myPos))])
        return "\t".join([self.name,self.switchStrand[myStrand],myChr,str(myPos[0]),'255',cigar,'*','0','0',self.read,self.qual,'NM:i:'+str(self.subs)])


    def jxnSearch(self):
        self.spliced = False
        myChr,myStrand,chrPos = self.hgLoc
        myGene,genePos        = self.geneLoc

        chrStrings  = [str(s) for s in chrPos]
        geneStrings = [str(s) for s in genePos]

        if len(chrPos) != len(genePos):
            print "WRF OMG NOT SAME LENGTH??? FUUUK"
            sys.exit()
        
        myOffset = 0 ; minLen = 10; mySplices = []; myPieces = []

        for i in range(0,len(chrPos),2):

            myDist = (chrPos[i+1] - chrPos[i]) + 1
            if myDist > minLen or (i>0 and i + 2 < len(chrPos)):

                geneData = myGene + ":" + ",".join(geneStrings[i:i+2])
                
                if i > 1:
                    mySplices.append("".join([myChr,":",chrStrings[i-1],"-",chrStrings[i],"|",myGene,':',geneStrings[i-1],'-',geneStrings[i]]))

                myPieces.append([self.name,myChr,myStrand,chrStrings[i],chrStrings[i+1],self.read[myOffset:myOffset+myDist],self.ref[myOffset:myOffset+myDist],self.qual,str(self.subs),geneData])
            myOffset += myDist 
        if len(mySplices) > 0:
            self.spliced = True
            totalSplice = ["EXONIC-SPLICED|"+",".join(mySplices)]
            for i in range(len(myPieces)):  myPieces[i] = myPieces[i] + totalSplice
        else:
            self.spliced = False
            for i in range(len(myPieces)):  myPieces[i] = myPieces[i] + ["EXONIC-CONT|"+myChr+"|"+myGene]
        self.spliceJxns = mySplices
        self.splitReads = myPieces


    def hgJxnSearch(self):

        myOffset = 0; minLen = 0; mySplice = None; myPieces = []
        
        myChr,myStrand,chrPos = self.hgLoc
        chrStrings  = [str(s) for s in chrPos]
        
        for i in range(0,len(chrPos),2):
            myDist = (chrPos[i+1] - chrPos[i]) + 1
            if myDist > minLen or (i>0 and i + 2 < len(chrPos)):
                if i > 1:   mySplice="".join(["INTERGENIC-SPLICED|",myChr,":",chrStrings[i-1],"-",chrStrings[i]])
                myPieces.append([self.name,myChr,myStrand,chrStrings[i],chrStrings[i+1],self.read[myOffset:myOffset+myDist],self.ref[myOffset:myOffset+myDist],self.qual,str(self.subs),"INTERGENIC"])
            myOffset += myDist 
        if mySplice:
            self.spliced = True
            for i in range(len(myPieces)):  myPieces[i].append(mySplice)
        else:
            self.spliced = False
            for i in range(len(myPieces)):  myPieces[i].append("INTERGENIC-CONT|"+myChr)
        self.spliceJxns = [mySplice]
        self.splitReads = myPieces

















































    def removeDuplicates(self):
        n=1; self.readLines.sort()
        while n < len(self.readLines):
            if self.readLines[n] == self.readLines[n-1]:
                self.readLines.remove(self.readLines[n])
            else:
                n+=1
        self.readLines.sort()

        
        if len(self.readLines) == 1:
            self.hgUniq = True; self.geneUniq = True
        elif self.readLines[0][0] == self.readLines[-1][0]:
            self.hgUniq  =  True
            self.geneUniq = False
        else:
            self.hgUniq   = False
            self.geneUniq = False


    def spliceEval(self):

        tmpGeneNames = set([])
        for i in range(len(self.readLines)):
            hgTmp = self.readLines[i][0]; geneTmp = self.readLines[i][1]; tmpGeneNames.add(geneTmp[0])
            if len(geneTmp[1]) > 2:
                geneHead = geneTmp[1][1] - geneTmp[1][0]
                geneTail = geneTmp[1][len(geneTmp[1])-1] - geneTmp[1][len(geneTmp[1])-2]
                read,ref = self.readLines[i][3][0],self.readLines[i][3][1]
                strands   = self.readLines[i][2]
                if geneHead < self.minOverlap:
                    geneTmp[1] = geneTmp[1][2::]
                    hgTmp[1]   = hgTmp[1][2::]
                    read = read[geneHead+1::]
                    ref  = ref[geneHead+1::]
                    
                if geneTail < self.minOverlap:
                    geneTmp[1] = geneTmp[1][0:len(geneTmp[1])-2]
                    hgTmp[1]   = hgTmp[1][0:len(hgTmp[1])-2]
                    read = read[0:len(read)-(geneTail+1)]
                    ref  = ref[0:len(ref)-(geneTail+1)]
                
                self.readLines[i] = [hgTmp,geneTmp,strands,[read,ref]]
        self.removeDuplicates()

        if not self.geneUniq and len(tmpGeneNames)==1:
            if self.readLines[0][0][1][0] == self.readLines[-1][0][1][0] and self.readLines[0][0][1][-1] == self.readLines[-1][0][1][-1]:
                ### CASE FALSE JUNCTION -- INTRONING JXN NEARBY --- ###
                if len(self.readLines[0][0][1]) < len(self.readLines[-1][0][1]):
                    MASTERmap = self.readLines[0]
                else:
                    MASTERmap = self.readLines[-1]
                for i in range(len(self.readLines)):
                    if self.readLines[i] != MASTERmap and self.readLines[i][0][1][0] == MASTERmap[0][1][0] and self.readLines[i][0][1][-1] == MASTERmap[0][1][-1]:
                        self.readLines[i] = MASTERmap

            self.removeDuplicates()
        


        



    def removePalindromes(self):


        if  self.readLines[0][0][1] != self.readLines[-1][0][1] and self.readLines[0][0][1] == self.readLines[-1][0][1][-1::-1]:

            tHg = self.readLines[0][0]; tmpSeqs = self.readLines[0][3]; tLocs = self.readLines[0][0][1]
            for i in range(len(self.readLines)):
                if self.readLines[i][0][1] == tLocs[-1::-1]:
                    if self.readLines[i][2][0] == "+":
                        newStr = '-'
                    else:
                        newStr = '+'
                    self.readLines[i] = [tHg,self.readLines[i][1],(newStr,self.readLines[i][2][1]),tmpSeqs]


    def evaluate(self):

        self.removeDuplicates()

        if self.strandRule != None:
            self.removeAntiSense()
        else:
            self.sense = True

        self.spliceEval()

        if not self.hgUniq:
            self.removePalindromes()
            self.removeDuplicates()
        hgLocs = set([]); geneLocs = set([]); seqs = set([]); strands = set([]); features = set([])
        for i in range(len(self.readLines)):
            hgLocs.add((self.readLines[i][0][0],self.readLines[i][2],tuple(self.readLines[i][0][1])))
            geneLocs.add((self.readLines[i][1][0],self.readLines[i][2],tuple(self.readLines[i][1][1])))
            seqs.add(tuple(self.readLines[i][3]))
            strands.add(self.readLines[i][2][0])

        self.hgLocs = list(hgLocs)
        self.geneLocs = list(geneLocs)
        self.seqs = list(seqs)
        mapStrands = list(strands)


        if self.fType == "INTERGENIC":

            if len(hgLocs) > 1:
                self.valid = False
            else:
                self.seq = self.seqs[0][0]; self.ref = self.seqs[0][1]

        elif len(self.seqs) > 1 or len(mapStrands)>1:

            self.valid = False
            multiGenome = set([])
            for i in range(len(self.readLines)):
                multiGenome.add((self.readLines[i][0][0],self.readLines[i][2],tuple(self.readLines[i][0][1]),tuple(self.readLines[i][3])))

            self.multiGenome = list(multiGenome)
        else:
            
            ## NOTE FIX THE QUAL REVERSE IT WHEN NECESSARY ##
            self.seq = self.seqs[0][0]; self.ref = self.seqs[0][1]
            if mapStrands.pop() == "-":
                self.qual = self.qual[-1::-1]
            if len(self.qual) > len(self.seq):
                self.qual = self.seq[0:len(self.seq)]
            
########################################################################################


    def findSplices(self):
        self.spliceTags = []
        if self.sense and self.geneUniq and self.valid and self.hgUniq:
            if len(self.geneLocs[0][2]) > 2:
                self.spliced = True
                tmpGene = self.geneLocs[0][2]; tmpHg = self.hgLocs[0][2];
                geneJumps = []; hgJumps = []; tmpName = self.geneLocs[0][0]; tmpChr = self.hgLocs[0][0]; tmpStrand = self.hgLocs[0][1][1]
                for i in range(1,len(tmpGene),2):
                    if i+1 < len(tmpGene):
                        geneJumps.append((str(tmpGene[i]),str(tmpGene[i+1])))
                        hgJumps.append((str(tmpHg[i]),str(tmpHg[i+1])))
                for j in range(len(geneJumps)):
                    junctionID = tmpChr+"_"+tmpStrand+"_"+hgJumps[j][0]+","+hgJumps[j][1]+"|"+tmpName+"_"+geneJumps[j][0]+","+geneJumps[j][1]
                    self.spliceTags.append(junctionID)
        if len(self.spliceTags)==0:
            self.spliceTags.append("UNSPLICED")


    def findGeneMaps(self):
        
        if not self.valid or self.fType == "INTERGENIC":
            self.geneMaps =  []
        elif self.geneUniq:
            self.geneMaps = [self.geneLocs[0][0]]
        else:
            self.geneMaps = tuple(sorted(list(set([self.geneLocs[i][0] for i in range(len(self.geneLocs))]))))

    def roundChrMaps(self):
        roundMaps=[]        
        tmpChr = self.hgLocs[0][0]
        tmpLocs = self.hgLocs[0][2]
        for i in range(0,len(tmpLocs),2):
            roundMaps.append(tmpLocs[i]-tmpLocs[i]%1000)
        return tmpChr,roundMaps


    def splitLocations(self):

        if self.valid and self.hgUniq and self.sense:
            hg = self.hgLocs[0]
                 
            if len(hg[2])>2: self.spliced = True
            else:            self.spliced = False
           

            if self.fType == "INTRONIC":
                if self.spliced: SPLICE_INFO="INTRONIC-JXNCROSS|"+hg[0]+":"
                else:       SPLICE_INFO="INTRONIC-CONTAINED|"+hg[0]+":"
            elif self.fType == "INTERGENIC":
                if self.spliced: SPLICE_INFO="INTERGENIC-SPLICED|"+hg[0]+":"
                else:       SPLICE_INFO="INTERGENIC-UNSPLICED|"+hg[0]+":"
            else:
                if self.spliced: SPLICE_INFO="EXONIC-SPLICED|"+hg[0]+":"
                else:       SPLICE_INFO="EXONIC-UNSPLICED|"+hg[0]+":"

            ## ------------------------------------------------------------------------------------------------------------------- ##

            hgTmp, tmpSplices,spliceData = [],[],[]
            tmpStart,n = 0,0
              
            for k in range(0,len(hg[2]),2):
                hgStart = hg[2][k]; hgEnd = hg[2][k+1]
                hgTmp.append([hgStart,hgEnd,self.seq[tmpStart:tmpStart+(hgEnd-hgStart)+1],self.ref[tmpStart:tmpStart+(hgEnd-hgStart)+1],self.qual[tmpStart:tmpStart+(hgEnd-hgStart)+1]])
                tmpStart += (hgEnd-hgStart) +1
                if k>1:     tmpSplices.append(str(hg[2][k-1])+"-"+str(hg[2][k]))
            if len(tmpSplices) !=0:
                SPLICE_INFO+=",".join(tmpSplices)
            else:
                SPLICE_INFO+="CONT"

            geTmp = [[] for t in hgTmp]
            for ge in self.geneLocs:
                tmpSplices = []
                for k in range(0,len(ge[2]),2):
                    geTmp[n].append(ge[0]+":"+str(ge[2][k])+","+str(ge[2][k+1]))
                    if k>1: tmpSplices.append((str(ge[2][k-1])+"-"+str(ge[2][k])))
                    n+=1
                if len(tmpSplices) != 0:
                    spliceData.append(SPLICE_INFO+"|"+ge[0]+":"+",".join(tmpSplices))
                else:
                    spliceData.append(SPLICE_INFO+"|"+ge[0]+":CONT")
                n=0

#            for i in range(len(hgTmp)):
#                print self.ID,hg[0],hg[1][0]," ".join([str(s) for s in hgTmp[i]]),self.subs,"|"," ".join(geTmp[i]),"|"," ".join(spliceData)






































































                

                    



    def writeSamHeader(self,outfile):
        
        hgHeader=["@HD\tVN:0.1.5c\tSO:queryname",
                "@SQ\tSN:chr10\tLN:135534748",
                "@SQ\tSN:chr11\tLN:135006517",
                "@SQ\tSN:chr12\tLN:133851896",
                "@SQ\tSN:chr13\tLN:115169879",
                "@SQ\tSN:chr14\tLN:107349541",
                "@SQ\tSN:chr15\tLN:102531393",
                "@SQ\tSN:chr16\tLN:90354754",
                "@SQ\tSN:chr17\tLN:81195211",
                "@SQ\tSN:chr18\tLN:78077249",
                "@SQ\tSN:chr19\tLN:59128984",
                "@SQ\tSN:chr1\tLN:249250622",
                "@SQ\tSN:chr20\tLN:63025521",
                "@SQ\tSN:chr21\tLN:48129896",
                "@SQ\tSN:chr22\tLN:51304567",
                "@SQ\tSN:chr2\tLN:243199374",
                "@SQ\tSN:chr3\tLN:198022431",
                "@SQ\tSN:chr4\tLN:191154277",
                "@SQ\tSN:chr5\tLN:180915261",
                "@SQ\tSN:chr6\tLN:171115068",
                "@SQ\tSN:chr7\tLN:159138664",
                "@SQ\tSN:chr8\tLN:146364023",
                "@SQ\tSN:chr9\tLN:141213432",
                "@SQ\tSN:chr17_ctg5_hap1\tLN:1680829",
                "@SQ\tSN:chr4_ctg9_hap1\tLN:590427",
                "@SQ\tSN:chr6_apd_hap1\tLN:4622291",
                "@SQ\tSN:chr6_cox_hap2\tLN:4795372",
                "@SQ\tSN:chr6_dbb_hap3\tLN:4610397",
                "@SQ\tSN:chr6_mann_hap4\tLN:4683264",
                "@SQ\tSN:chr6_mcf_hap5\tLN:4833399",
                "@SQ\tSN:chr6_qbl_hap6\tLN:4611985",
                "@SQ\tSN:chr6_ssto_hap7\tLN:4928568",
                "@SQ\tSN:chrM\tLN:16572",
                "@SQ\tSN:chr11_gl000202_random\tLN:40104",
                "@SQ\tSN:chr17_gl000203_random\tLN:37499",
                "@SQ\tSN:chr17_gl000204_random\tLN:81311",
                "@SQ\tSN:chr17_gl000205_random\tLN:174589",
                "@SQ\tSN:chr17_gl000206_random\tLN:41002",
                "@SQ\tSN:chr18_gl000207_random\tLN:4263",
                "@SQ\tSN:chr19_gl000208_random\tLN:92690",
                "@SQ\tSN:chr19_gl000209_random\tLN:159170",
                "@SQ\tSN:chr1_gl000191_random\tLN:106434",
                "@SQ\tSN:chr1_gl000192_random\tLN:547497",
                "@SQ\tSN:chr21_gl000210_random\tLN:27683",
                "@SQ\tSN:chr4_gl000193_random\tLN:189790",
                "@SQ\tSN:chr4_gl000194_random\tLN:191470",
                "@SQ\tSN:chr7_gl000195_random\tLN:182897",
                "@SQ\tSN:chr8_gl000196_random\tLN:38915",
                "@SQ\tSN:chr8_gl000197_random\tLN:37176",
                "@SQ\tSN:chr9_gl000198_random\tLN:90086",
                "@SQ\tSN:chr9_gl000199_random\tLN:169875",
                "@SQ\tSN:chr9_gl000200_random\tLN:187036",
                "@SQ\tSN:chr9_gl000201_random\tLN:36149",
                "@SQ\tSN:chrUn_gl000211\tLN:166567",
                "@SQ\tSN:chrUn_gl000212\tLN:186859",
                "@SQ\tSN:chrUn_gl000213\tLN:164240",
                "@SQ\tSN:chrUn_gl000214\tLN:137719",
                "@SQ\tSN:chrUn_gl000215\tLN:172546",
                "@SQ\tSN:chrUn_gl000216\tLN:172295",
                "@SQ\tSN:chrUn_gl000217\tLN:172150",
                "@SQ\tSN:chrUn_gl000218\tLN:161148",
                "@SQ\tSN:chrUn_gl000219\tLN:179199",
                "@SQ\tSN:chrUn_gl000220\tLN:161803",
                "@SQ\tSN:chrUn_gl000221\tLN:155398",
                "@SQ\tSN:chrUn_gl000222\tLN:186862",
                "@SQ\tSN:chrUn_gl000223\tLN:180456",
                "@SQ\tSN:chrUn_gl000224\tLN:179694",
                "@SQ\tSN:chrUn_gl000225\tLN:211174",
                "@SQ\tSN:chrUn_gl000226\tLN:15009",
                "@SQ\tSN:chrUn_gl000227\tLN:128375",
                "@SQ\tSN:chrUn_gl000228\tLN:129121",
                "@SQ\tSN:chrUn_gl000229\tLN:19914",
                "@SQ\tSN:chrUn_gl000230\tLN:43692",
                "@SQ\tSN:chrUn_gl000231\tLN:27387",
                "@SQ\tSN:chrUn_gl000232\tLN:40653",
                "@SQ\tSN:chrUn_gl000233\tLN:45942",
                "@SQ\tSN:chrUn_gl000234\tLN:40532",
                "@SQ\tSN:chrUn_gl000235\tLN:34475",
                "@SQ\tSN:chrUn_gl000236\tLN:41935",
                "@SQ\tSN:chrUn_gl000237\tLN:45868",
                "@SQ\tSN:chrUn_gl000238\tLN:39940",
                "@SQ\tSN:chrUn_gl000239\tLN:33825",
                "@SQ\tSN:chrUn_gl000240\tLN:41934",
                "@SQ\tSN:chrUn_gl000241\tLN:42153",
                "@SQ\tSN:chrUn_gl000242\tLN:43524",
                "@SQ\tSN:chrUn_gl000243\tLN:43342",
                "@SQ\tSN:chrUn_gl000244\tLN:39930",
                "@SQ\tSN:chrUn_gl000245\tLN:36652",
                "@SQ\tSN:chrUn_gl000246\tLN:38155",
                "@SQ\tSN:chrUn_gl000247\tLN:36423",
                "@SQ\tSN:chrUn_gl000248\tLN:39787",
                "@SQ\tSN:chrUn_gl000249\tLN:38503",
                "@SQ\tSN:chrX\tLN:155270561",
                "@SQ\tSN:chrY\tLN:59373567"
                "@RG\tID:knowles.fastq\tSM:knowles.fastq",
                "@PG\tID:PerM\tVN:0.4.0"]

        for h in hgHeader: outfile.write("%s\n" % h)



    



##############################################################################################################





