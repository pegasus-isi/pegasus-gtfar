#!/usr/bin/env python


import sys
from MapRead import *
from Sequence import *
from Utilities import *



from collections import defaultdict as dd

#from ..gtTools.seq import *
#from ..gtTools.util import *
from math import fabs

##########################################################################################################################################
#####################################################  MAPFILE-FILE CLASS START ##########################################################
##########################################################################################################################################

class MapFile:

    def __init__(self,fileHandle,prefix,FTYPE,strandSpecific=None):

        ## FIGURE OUT THE FILE FORMAT ##
        self.fname = open(fileHandle); self.prefix = prefix; self.fType = FTYPE; self.sense = True; self.showAMBIG = True; self.index =0; self.open = True; self.antiSense = None; self.key = {}; self.FULLSEQ=False

        if fileHandle.split(".")[-1] == "mapping":      self.format = "MAPPING"
        elif fileHandle.split(".")[-1] == "sam":        self.format = "SAM"
        else:                                           errorQuit('Invalid mapfile')

        self.firstLine = self.fname.readline().split()
        if self.firstLine==[]: errorQuit('Empty mapfile')
           
        if self.format == "MAPPING":
            if strandSpecific != None:
                if strandSpecific == '0' or strandSpecific == '+':      self.antiSense = '-'
                else:                                                   self.antiSense = '+'
        else:
            if strandSpecific!= None:
                if strandSpecific == '0' or strandSpecific == "+":      self.antiSense = '16'
                else:                                                   self.antiSense = '0'
        
        self.samOut, self.mapOut,self.AntiOut, self.geneOut,self.spliceOut,self.featOut = None, None , None, None, None, None

##############################################################################################################
############################################   ADD KEY   #####################################################
##############################################################################################################

    def loadKey(self,keyFile):    
        if self.fType == "EXONIC":
            for line in open(keyFile):
                line=line.split()
                refType=line[0].split("|")[6][0:4]
                if refType=="GENE" or refType=="TRAN":
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]
        elif self.fType == "INTRONIC":
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
        elif self.fType == "GAPPED":
            for line in open(keyFile):
                line = line.split()
                refType=line[0].split("|")[6][0:4]
                if refType=="GENE":
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]
        elif self.fType == "HG19":
            return
        else: 
            errorQuit("Invalid Ref Type")





















##############################################################################################################
############################################   READ RELOCATION  ##############################################
##############################################################################################################


    @staticmethod
    def relocate(fPos,key,DIST):

        fKey,gKey,hKey = key[0],key[1],key[2];  POST = False; gPos =[]; hPos =[]
        
        for i in range(len(fKey)):
            
            if POST or (fPos >= fKey[i][0] and fPos <= fKey[i][1]):

                if not POST:
                    myOffset = fPos - fKey[i][0]
                gPos.extend([gKey[i][0]+ myOffset, gKey[i][0] + myOffset + DIST  ])

                if hKey[i][0] < hKey[i][1]:
                   hPos.extend( [ hKey[i][0] + myOffset, hKey[i][0] + myOffset + DIST ] )
                else:
                   hPos.extend( [ hKey[i][0] - myOffset, hKey[i][0] - myOffset - DIST ] )

                if fKey[i][1] >= fPos + DIST and fKey[i][1] >= fKey[i][0] + DIST:
                    break
                else:
                    gPos[-1] = gKey[i][1]; hPos[-1] = hKey[i][1]
                    DIST -= (1 + fKey[i][1] - (fKey[i][0] + myOffset) )
                    POST = True; myOffset =0
        return gPos,hPos


    @staticmethod
    def relocateGenePos(fPos,key,DIST):
        fKey,gKey,hKey = key[0],key[1],key[2];  POST = False; gPos =[]; hPos =[]; eData = []
        ### DO THIS LIKE A SIMPLETON ###


        for i in range(len(gKey)):

            if POST or (fPos>=gKey[i][0] and fPos<= gKey[i][1]):

                if not POST:  myOffset = fPos - gKey[i][0]

                if fPos == gKey[i][0]: eData.append("EX-STR")
                else:                  eData.append("EX-MID")

                gPos.extend([gKey[i][0] + myOffset, gKey[i][0] + myOffset + DIST])

                if hKey[i][0] < hKey[i][1]: hPos.extend( [ hKey[i][0] + myOffset, hKey[i][0] + myOffset + DIST ] )
                else:                       hPos.extend( [ hKey[i][0] - myOffset,  hKey[i][0] - myOffset - DIST ] )

                if gKey[i][1] >= fPos + DIST and gKey[i][1] >= gKey[i][0] + DIST:

                    if fPos + DIST == gKey[i][1]:
                        eData.append("EX-END")
                    else:
                        eData.append("EX-MID")
                    break
                else:
                    eData.append("EX-END")
                    gPos[-1] = gKey[i][1]; hPos[-1] = hKey[i][1]
                    DIST -= (1 + gKey[i][1] - (gKey[i][0] + myOffset) )
                    POST = True; myOffset = 0                 

            elif  (i>0)  and  (fPos>gKey[i-1][1] and fPos < gKey[i][0]) :



                hgPrev = hKey[i-1]; myOffset = fPos - gKey[i-1][1];    gPos.extend([ fPos, fPos+DIST ])

                eData.append("INT")

                if hgPrev[0] < hgPrev[1]: hPos.extend( [hgPrev[1] + myOffset, hgPrev[1] + myOffset + DIST ] )
                else:                     hPos.extend( [hgPrev[1] - myOffset, hgPrev[1] - myOffset - DIST ] )

                if fPos + DIST < gKey[i][0] and gKey[i][0] >= gKey[i-1][1] + DIST:
                    if fPos + DIST == gKey[i][0] -1: eData.append("IN-END")
                    else:                            eData.append("IN-MID")
                    break
                else:
                    eData.extend(["IN-END","EX-STR"])
                    gPos[-1] = gKey[i][0]; hPos[-1] = hKey[i][0]
                    DIST -= (1 + gKey[i][0] - (gKey[i-1][1] + myOffset))
                    POST = True; myOffset = 0

        return gPos,hPos,eData


    def relocateCigar(self,fPos,cigar,key,FULLSEQ=False):
        fPos = fPos-1
        cigar = cigar.strip().split('M')
        
        if not FULLSEQ:
            gInit, hInit  = self.relocate(fPos,key,int(cigar[0])-1)
            nPos = fPos + int(cigar[0])

            for pair in cigar[1:len(cigar)-1]:
                c=pair.split('N')
                gap = int(c[0]); match = int(c[1])
                gNext,hNext = self.relocate(nPos + gap, key, match-1) 
                gInit.extend(gNext); hInit.extend(hNext)
                nPos = nPos + gap + match  ## COULD BE WRONG ##
            return gInit,hInit
        else:


            gInit, hInit , eData = self.relocateGenePos(fPos,key,int(cigar[0])-1)
            nPos = fPos + int(cigar[0])

            for pair in cigar[1:len(cigar)-1]:
                c=pair.split('N')
                gap = int(c[0]); match = int(c[1])
                gNext,hNext, eNext = self.relocateGenePos(nPos + gap, key, match-1) 
                gInit.extend(gNext); hInit.extend(hNext); eData.append("GAP"), eData.extend(eNext)
                nPos = nPos + gap + match  ## COULD BE WRONG ##
            return gInit,hInit





##############################################################################################################
############################################   PARSING READ   ################################################
##############################################################################################################


    def parseExonicLines(self,mapLines,strandCnt,sense):

            tMaps =[]; tCats = [];  catLines = []; newLines = []; sense = True

            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:

                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]
                    geneTmp,hgTmp = self.relocate(int(mapLines[i][3]),tmpKey,len(mapLines[i][1])-1)
                    hgStrand,hgLoc = trueStrand(mapLines[i][5],s[2], hgTmp )


             #       print "OK",geneTmp,hgTmp
             #       print "FUCK",hgStrand,hgLoc


                    myLine = [ ( s[3] , hgStrand, hgLoc  ) , ( s[0],mapLines[i][5],geneTmp ) , ( mapLines[i][1],mapLines[i][4] ) ]

                    if s[len(s)-1]=="CATSEQ":
                        tCats.append(mapLines[i][2])
                        if catLines == []:
                            ID = mapLines[i][0]; qual = mapLines[i][8]; subs= int(mapLines[i][6]); catLines.append(myLine)
                        else:
                            if myLine not in catLines:      catLines.append(myLine)
                    
                    else:
                        tMaps.append(mapLines[i][2])
                        if newLines == []:
                     #       print "PASS",mapLines[i][2],geneTmp,hgTmp
                            ID = mapLines[i][0]; qual = mapLines[i][8]; subs= int(mapLines[i][6]); newLines.append(myLine)
                        else:
                            if myLine not in newLines:
                      #          print "PASS",mapLines[i][2],geneTmp,hgTmp
                                newLines.append(myLine)

            if len(newLines) > 0: self.read = MapRead(self.format,sense,ID,qual,subs,tMaps,newLines)
            else:                 self.read = MapRead(self.format,sense,ID,qual,subs,tCats,catLines)
            

##############################################################################################################

    def parseIntronicLines(self,mapLines,strandCnt,sense):

            tMaps =[]; tCats = [];  catLines = []; newLines = []; sense = True

            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:

                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]
                    geneTmp,hgTmp = self.relocate(int(mapLines[i][3]),tmpKey,len(mapLines[i][1])-1)
                    


                    hgStrand,hgLoc = trueStrand(mapLines[i][5],s[2], hgTmp )

                    myLine = [ ( s[3] , hgStrand, hgLoc  ) , ( s[0],mapLines[i][5],geneTmp ) , ( mapLines[i][1],mapLines[i][4] ) ]

                    if s[len(s)-1][0:4]=="FLNK":
                        tCats.append(mapLines[i][2])
                        if catLines == []:
                            ID = mapLines[i][0]; qual = mapLines[i][8]; subs= int(mapLines[i][6]); catLines.append(myLine)
                        else:
                            if myLine not in catLines:      catLines.append(myLine)
                    
                    else:
                        tMaps.append(mapLines[i][2])
                        if newLines == []:
                            ID = mapLines[i][0]; qual = mapLines[i][8]; subs= int(mapLines[i][6]); newLines.append(myLine)
                        else:
                            if myLine not in newLines:
                                newLines.append(myLine)
             
            if len(newLines) > 0: self.read = MapRead(self.format,sense,ID,qual,subs,tMaps,newLines)
            else:                 self.read = MapRead(self.format,sense,ID,qual,subs,tCats,catLines)
            

##############################################################################################################

    def parseGappedLines(self,mapLines,strandCnt,sense):


            ### NOTICE WE ARE TACILITY ASSUMING SAM FORMAT ###


            tMaps =[]; tCats = [];  catLines = []; newLines = []; sense = True

            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:

                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]

                    if s[-1] == "FULLSEQ":
                        self.FULLSEQ = True
                        geneTmp, hgTmp = self.relocateCigar(int(mapLines[i][3]),mapLines[i][5],tmpKey,FULLSEQ=True)
                    else:
                        geneTmp, hgTmp = self.relocateCigar(int(mapLines[i][3]),mapLines[i][5],tmpKey,)

                    
                    signStrand = sam2Map(mapLines[i][1])
                    

                    hgStrand,hgLoc = trueStrand(signStrand,s[2], hgTmp)
                    myLine = [ ( s[3] , hgStrand, hgLoc  ) , ( s[0],signStrand,geneTmp ) , ( mapLines[i][9],mapLines[i][9] ) ]

                    if s[len(s)-1]=="CATSEQ":
                        tCats.append(mapLines[i][2])
                        if catLines == []:
                            ID = mapLines[i][0]; qual = mapLines[i][10]; subs= int(mapLines[i][11].split(':')[-1]); catLines.append(myLine)
                        else:
                            if myLine not in catLines:      catLines.append(myLine)
                    else:
                        tMaps.append(mapLines[i][2])
                        if newLines == []:
                            ID = mapLines[i][0]; qual = mapLines[i][10]; subs= int(mapLines[i][11].split(':')[-1]); newLines.append(myLine)
                        else:
                            if myLine not in newLines:
                                newLines.append(myLine)

            if len(newLines) > 0: self.read = MapRead(self.format,sense,ID,qual,subs,tMaps,newLines)
            else:                 self.read = MapRead(self.format,sense,ID,qual,subs,tCats,catLines)
            

##############################################################################################################


    def parseGenomicLines(self,mapLines,strandCnt,sense):
        
        if self.index == 1:
            self.spliceExp = dd(int)
        if len(mapLines) == 1:
            hgData = (mapLines[0][2], sam2Map(mapLines[0][1]), cigar2List(int(mapLines[0][3]),mapLines[0][5])) 
            
            spliceStr=hgData[0]+":"
            for i in range(1,len(hgData[2])-1,2):
                spliceStr+=str(hgData[2][i])+">"+str(hgData[2][i+1])

            self.spliceExp[spliceStr]+=1


        

##############################################################################################################




























##############################################################################################################
############################################   GETTING READ   ################################################
##############################################################################################################

    def getReads(self):

        firstID = self.firstLine[0]; strandSet = set([]); mapLines=[]; self.index+=1
        if self.format == "MAPPING":
            while self.firstLine[0] == firstID:
                mapLines.append(self.firstLine)
                strandSet.add(self.firstLine[5])
                self.firstLine = self.fname.readline().split()
                if self.firstLine == []:
                    self.open = False
                    break

            strandCnt = len(strandSet); sense = True
            
            if strandCnt == 1:
                if strandSet.pop() == self.antiSense: sense = False

            if self.fType == "EXONIC":
                self.parseExonicLines(mapLines,strandCnt,sense)
            
            elif self.fType == "INTRONIC":
                self.parseIntronicLines(mapLines,strandCnt,sense)

        else:
            while self.firstLine[0] == firstID:
                mapLines.append(self.firstLine)
                strandSet.add(self.firstLine[1])
                self.firstLine = self.fname.readline().split()
                if self.firstLine == []:
                    self.open = False; break
            strandCnt = len(strandSet); sense = True
           
            if strandCnt == 1:
                if strandSet.pop() == self.antiSense: sense = False
            
            if self.fType == "GAPPED":
                self.parseGappedLines(mapLines,strandCnt,sense)
            elif self.fType == "HG19":
                self.parseGenomicLines(mapLines,strandCnt,sense)














































##############################################################################################################
############################################   STORING DATA   ################################################
##############################################################################################################


    @staticmethod
    def findSpliceJxns(hg,gene):
        minOverlap=5; geneJxns=[]; hgJxns =[]
        if len(gene[2])<4:
            return None
        else:
            for i in range(0,len(gene[2])-2,2):
                if gene[2][i+1]-gene[2][i] > minOverlap and gene[2][i+3]-gene[2][i+2] > minOverlap:
                    geneJxns.append(str(gene[2][i+1])+">"+str(gene[2][i+2]))
                    hgJxns.append(str(hg[2][i+1])+">"+str(hg[2][i+2]))
            
            if len(hgJxns)==0:
                return None
            else:
                spliceStr=hg[0]+':'+','.join(hgJxns)+'|'+gene[0]+'='+'.'.join(geneJxns)
                return spliceStr

#############################################################################################################################################################


    def storeExpression(self):

        if self.index == 1:
            self.geneExp = {}
            self.featExp = dd(int)
            self.spliceExp = dd(int); self.geneSeen = set([])
       

        for i in range(len(self.read.hgLocs)):
            for gene in self.read.geneLocs[i]:
                if gene[0] not in self.geneExp:
                    self.geneExp[gene[0]]= [0,0,0]; self.geneSeen.add(gene[0])
                
                if self.read.uniq and not self.read.multi:
                    self.geneExp[gene[0]][0]+=1
                    jxn = self.findSpliceJxns(self.read.hgLocs[0],gene)
                    if jxn:     self.spliceExp[jxn]+=1

                    for f in self.read.features:
                        self.featExp[f]+=1
                elif self.read.uniq:
                    self.geneExp[gene[0]][1]+=1
                else:
                    self.geneExp[gene[0]][2]+=1























##############################################################################################################
####################################   PRINTING DATA #########################################################
##############################################################################################################
                
    def writeNovelGenes(self):

        self.spliceOut = open(self.prefix+'_splice.cnts','w')
        for s in self.spliceExp.keys():
            self.spliceOut.write("%s NOVEL %s\n" % (s,self.spliceExp[s]))

    def writeExpression(self):

        self.geneOut = open(self.prefix+'_gene.cnts','w')
        self.spliceOut = open(self.prefix+'_splice.cnts','w')
        self.featOut = open(self.prefix+'_feat.cnts','w')

        if self.fType == "GAPPED":
            if self.FULLSEQ: self.fType = "GENE-GAP"
            else:            self.fType = "CAT-GAP" 

        for key in self.key.keys():
            key = key.split("|")
            if key[6] == "GENE":
                k="|".join(key[0:6])
                exSize = int(key[7]); gSize = int(key[8])
                if key[0] in self.geneSeen:
                    self.geneOut.write("%s %s %s | %s %s %s %s\n" % (k,exSize,gSize,self.fType,self.geneExp[key[0]][0],self.geneExp[key[0]][1],self.geneExp[key[0]][2]))
                else:
                    self.geneOut.write("%s %s %s | %s %s %s %s\n" % (k,exSize,gSize,self.fType,0,0,0))
        
        for s in self.spliceExp.keys():
            self.spliceOut.write("%s %s %s\n" % (s,self.fType,self.spliceExp[s]))

        for f in self.featExp.keys():
            self.featOut.write("%s %s %s\n" % (f,self.fType,self.featExp[f]))




    def writeLocations(self):

        self.read.visualizeStrings()

        if self.index == 1:
            self.samOut = open(self.prefix+'_vis.sam','w')
            self.mapOut = open(self.prefix+'_gene.loc','w')
            self.AntiOut = open(self.prefix+'_anti.loc','w')

        for s in self.read.samStrings:

            if self.showAMBIG or len(self.read.samStrings) == 1:
                self.samOut.write('%s\n' % s)

        for g in self.read.locStrings:
        
            if self.read.sense == False:

                self.AntiOut.write('%s\n' %g)
            else:
                if len(self.read.locStrings) == 1 or self.showAMBIG:
                
                    self.mapOut.write('%s\n' % g)


    def close(self):

        self.fname.close()
        if self.samOut: self.samOut.close()
        if self.mapOut: self.mapOut.close()
        if self.AntiOut: self.AntiOut.close()
        if self.geneOut: self.geneOut.close()
        if self.spliceOut: self.spliceOut.close()
        if self.featOut: self.featOut.close()








