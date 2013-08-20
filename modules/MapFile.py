#!/usr/bin/env python


import sys
from MapRead import *
from Sequence import *
from Utilities import *

from collections import defaultdict as dd
from math import fabs

##########################################################################################################################################
#####################################################  MAPFILE-FILE CLASS START ##########################################################
##########################################################################################################################################

class MapFile:

    def __init__(self,mapFile,keyFile,refType,GAPPED,prefix,strandSpecific=None):
 

        self.mapFile,self.keyFile,self.refType,self.GAPPED,self.prefix = mapFile,keyFile,refType,GAPPED,prefix

        self.mapHandle = open(self.mapFile); self.keyHandle = open(self.keyFile)
        
        self.invalid, self.showAMBIG, self.FULLSEQ, self.FLANKSEQ = False, False, False, False
        
        self.samOut, self.mapOut,self.AntiOut, self.geneOut,self.spliceOut,self.featOut,self.antiSense = None, None , None, None, None, None, None
        
        self.index = 0; self.key = {}

        self.sense, self.open = True, True

        ## CHECK: NON EMPTY FILE/ FILE FORMAT / STRAND-SPECIFICITY ##
        
        self.firstLine = self.mapHandle.readline().split()
        if self.firstLine==[]:  self.invalid = "EMPTY_MAPFILE"

        elif self.mapFile.split(".")[-1] == "mapping":
            self.format = "MAPPING"
            if strandSpecific != None:
                if strandSpecific == '0' or strandSpecific == '+':      self.antiSense = '-'
                else:                                                   self.antiSense = '+'


        elif self.mapFile.split(".")[-1] == "sam":        
            self.format = "SAM"
            if strandSpecific!= None:
                if strandSpecific == '0' or strandSpecific == "+":      self.antiSense = '16'
                else:                                                   self.antiSense = '0'
        else:
            self.invalid = "MAPFILE_EXTENSION"
        
        ##########################   LOAD THE KEY   #####################

        if self.invalid: return

 

        if self.refType   == "EXONIC":                               FEATURES = ["TRAN","GENE"]
        elif self.refType == "INTRONIC":                             FEATURES = ["ITRN","GENE","FLNK"]
        else:                                                        FEATURES = ["GENE"]
            
        for line in self.keyHandle:
            line    = line.split()
            feature = line[0].split("|")[6][0:4]
            if feature in FEATURES:
                if feature == "FLNK":
                    for i in range(1,len(line)):
                        line[i]=[[x.split("-")[0],x.split("-")[1]] for x in line[i].split(",")]
                        for j in range(len(line[i])):
                            for k in range(len(line[i][j])):
                                if line[i][j][k][0]=="N":       line[i][j][k]=int(line[i][j][k][1::])*-1
                                else:                           line[i][j][k]=int(line[i][j][k])
                    self.key[line[0]]=[line[1],line[2],line[3]]
                else:
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]

        
        





















        
        

##############################################################################################################
############################################   ADD KEY   #####################################################
##############################################################################################################



#print

    def loadKey(self,keyFile):    
        if self.refType == "EXONIC":
            for line in open(keyFile):
                line=line.split()
                refType=line[0].split("|")[6][0:4]
                if refType=="GENE" or refType=="TRAN":
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]
        elif self.refType == "INTRONIC":
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
        elif self.refType == "GAPPED":
            for line in open(keyFile):
                line = line.split()
                refType=line[0].split("|")[6][0:4]
                if refType=="GENE":
                    for i in range(1,len(line)):
                        line[i]=[[int(x.split("-")[0]),int(x.split("-")[1])] for x in line[i].split(",")]
                    self.key[line[0]]=[line[1],line[2],line[3]]
        elif self.refType == "HG19":
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
    def relocateFullGene(fPos,key,dList,fStrand):
        fKey,gKey,hKey = key[0],key[1],key[2]; hgList = []; k=0; i=0


        while len(hgList) < len(dList):
            if i == len(key[1]):
                for n in range(k,len(dList)):
                        if fStrand == "-":
                            hgList.append(key[2][i-1][1] - (dList[n] - key[1][i-1][1]) )
                        else:
                            hgList.append(key[2][i-1][1] + (dList[n] - key[1][i-1][1]) )
                break



            if dList[k] < key[1][i][1]:
                for j in range(k,len(dList)):
                    if dList[j] <= key[1][i][0]:
                        if fStrand == "-":
                            hgList.append( key[2][i][0] + (key[1][i][0] - dList[j]) )
                        else:
                            hgList.append( key[2][i][0] - (key[1][i][0] - dList[j]) )

                    elif dList[j] >= key[1][i][0] and dList[j] <= key[1][i][1]:
                   
                        if fStrand == "-":
                            hgList.append( key[2][i][0] - (dList[j] - key[1][i][0]) )
                        else:
                            hgList.append( key[2][i][0] + (dList[j] - key[1][i][0]) )
                    else:
                        k=j
                        break
            i+=1
        return dList,hgList



    def relocateCigar(self,fPos,cigar,key,fStrand,FULLSEQ=False):
        fPos = fPos-1
        if self.refType == "GENES":
            dList = cigar2List(fPos,cigar)
            geneLocs, hgLocs = self.relocateFullGene(fPos,key,dList,fStrand)
            return geneLocs, hgLocs
          
        if not FULLSEQ:
            cigar = cigar.strip().split('M')
            gInit, hInit  = self.relocate(fPos,key,int(cigar[0])-1)
            nPos = fPos + int(cigar[0])

            for pair in cigar[1:len(cigar)-1]:
                c=pair.split('N')
                gap = int(c[0]); match = int(c[1])
                gNext,hNext = self.relocate(nPos + gap, key, match-1) 
                gInit.extend(gNext); hInit.extend(hNext)
                nPos = nPos + gap + match  ## COULD BE WRONG ##
            return gInit,hInit
            
            




##############################################################################################################################################
########################################  PARSING READ (MAPPING FILE) ########################################################################
##############################################################################################################################################


#######################################  EXONIC UNGAPPED  ##########################################################

    def parseExonicLines(self,mapLines,strandCnt,sense):
            tMaps =[]; tCats = [];  catLines = []; newLines = []; sense = sense
            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:
                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]
                    geneTmp,hgTmp = self.relocate(int(mapLines[i][3]),tmpKey,len(mapLines[i][1])-1)
                    
                    hgStrand,hgLoc = trueStrand(mapLines[i][5],s[2], hgTmp )
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
                            ID = mapLines[i][0]; qual = mapLines[i][8]; subs= int(mapLines[i][6]); newLines.append(myLine)
                        else:
                            if myLine not in newLines:
                                newLines.append(myLine)

            if len(newLines) > 0: self.read = MapRead(self.format,sense,ID,qual,subs,tMaps,newLines)
            else:                 self.read = MapRead(self.format,sense,ID,qual,subs,tCats,catLines)
            

#########################################################################################################################################


#######################################  INTRONIC UNGAPPED  ##########################################################


    def parseIntronicLines(self,mapLines,strandCnt,sense):

            tMaps =[]; tCats = [];  catLines = []; newLines = []

            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:
                    

                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]
                    
                    
                    geneTmp,hgTmp = self.relocate(int(mapLines[i][3]),tmpKey,len(mapLines[i][1])-1)
                    hgStrand,hgLoc = trueStrand(mapLines[i][5],s[2], hgTmp )

                    #print hgStrand, hgLoc


                    myLine = [ ( s[3] , hgStrand, hgLoc  ) , ( s[0],mapLines[i][5],geneTmp ) , ( mapLines[i][1],mapLines[i][4] ) ]

                    if s[len(s)-2][0:4]=="FLNK":
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
            else:
                self.FLANKSEQ = True
                self.read = MapRead(self.format,sense,ID,qual,subs,tCats,catLines,seqType="FLANKSEQ")
    


#########################################################################################################################################


#######################################  INTRONIC UNGAPPED  #############################################################################



    def parseUngappedGenomicLines(self,mapLines,strandCnt,sense):
       
        if self.index == 1:
            self.novelAreas = {}
        if len(mapLines) == 1:

            mapChr = mapLines[0][2]
            mapPos = int(mapLines[0][3])
            mapRound = mapPos-mapPos%1000
            if mapChr not in self.novelAreas:
                self.novelAreas[mapChr]=dd(int)
                
            self.novelAreas[mapChr][mapRound]+=1
           
                
             



##############################################################################################################################################
########################################   SAM PARSERS BELOW   ###############################################################################
##############################################################################################################################################





    def parseGappedLines(self,mapLines,strandCnt,sense):

            ### NOTICE WE ARE TACILITY ASSUMING SAM FORMAT ###

            tMaps =[]; tCats = [];  catLines = []; newLines = []

            for i in range(len(mapLines)):

                if strandCnt == 1 or mapLines[i][5] != self.antiSense:

                    s = mapLines[i][2].split("|");  tmpKey = self.key["|".join(s[0:len(s)-1])+"|"]

                    if s[-1] == "FULLSEQ":
                        self.FULLSEQ = True
                        geneTmp, hgTmp = self.relocateCigar(int(mapLines[i][3]),mapLines[i][5],tmpKey,s[2],FULLSEQ=True)
                    else:
                        geneTmp, hgTmp = self.relocateCigar(int(mapLines[i][3]),mapLines[i][5],tmpKey,s[2],FULLSEQ=False)

                    
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


    def parseGappedGenomicLines(self,mapLines,strandCnt,sense):
        
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
                self.firstLine = self.mapHandle.readline().split()
                if self.firstLine == []:
                    self.open = False
                    break

            strandCnt = len(strandSet); sense = True
                        
            if strandCnt == 1:
                if strandSet.pop() == self.antiSense: sense = False

            if self.refType == "EXONIC":
                self.parseExonicLines(mapLines,strandCnt,sense)
            
            elif self.refType == "INTRONIC":
                self.parseIntronicLines(mapLines,strandCnt,sense)

            elif self.refType == "HG19":
                self.parseUngappedGenomicLines(mapLines,strandCnt,sense)

        else:
            while self.firstLine[0] == firstID:
                mapLines.append(self.firstLine)
                strandSet.add(self.firstLine[1])
                self.firstLine = self.mapHandle.readline().split()
                if self.firstLine == []:
                    self.open = False; break
            strandCnt = len(strandSet); sense = True
           
            if strandCnt == 1:
                if strandSet.pop() == self.antiSense: sense = False
             
            if self.GAPPED:
                
                if self.refType == "HG19":
                    self.parseGappedGenomicLines(mapLines,strandCnt,sense)
                if self.refType == "GENES" or self.refType == "EXONIC":

                    self.parseGappedLines(mapLines,strandCnt,sense)

            else:
                print self.refType
                errorQuit("SAM FORMAT FOR REQUIRES GAPPED ALIGNMENT")












































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


    def storeGeneCnts(self):

        if self.index == 1:
            self.geneExp = {}
            self.geneSeen = set([])
       

        for i in range(len(self.read.hgLocs)):
            for gene in self.read.geneLocs[i]:
                if gene[0] not in self.geneExp:
                    self.geneExp[gene[0]]= [0,0,0,0,0,0]; self.geneSeen.add(gene[0])
                
                if self.read.sense:
                    if self.read.uniq and not self.read.multi:
                        self.geneExp[gene[0]][0]+=1
                    elif self.read.uniq:
                        self.geneExp[gene[0]][1]+=1
                    else:
                        self.geneExp[gene[0]][2]+=1
                else:
                    if self.read.uniq and not self.read.multi:
                        self.geneExp[gene[0]][3]+=1
                    elif self.read.uniq:
                        self.geneExp[gene[0]][4]+=1
                    else:
                        self.geneExp[gene[0]][5]+=1

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

    

    def writeExpression(self):

        self.geneOut = open(self.prefix+'_gene.cnts','w')
        self.spliceOut = open(self.prefix+'_splice.cnts','w')

        if self.refType == "GAPPED":
            if self.FULLSEQ: self.refType = "GENE-GAP"
            else:            self.refType = "CAT-GAP" 

        for key in self.key.keys():
            key = key.split("|")
            if key[6] == "GENE":
                k="|".join(key[0:6])
                exSize = int(key[7]); gSize = int(key[8])
                if key[0] in self.geneSeen:
                    self.geneOut.write("%s %s %s | %s %s %s %s\n" % (k,exSize,gSize,self.refType,self.geneExp[key[0]][0],self.geneExp[key[0]][1],self.geneExp[key[0]][2]))
                else:
                    self.geneOut.write("%s %s %s | %s %s %s %s\n" % (k,exSize,gSize,self.refType,0,0,0))
        
        for s in self.spliceExp.keys():
            self.spliceOut.write("%s %s %s\n" % (s,self.refType,self.spliceExp[s]))

        if self.refType != "GAPPED":
            self.featOut = open(self.prefix+'_feat.cnts','w')

            for f in self.featExp.keys():
                self.featOut.write("%s %s %s\n" % (f,self.refType,self.featExp[f]))



    def writeGeneCnts(self):

        self.geneOut = open(self.prefix+'_gene.cnts','w')

        for key in self.key.keys():
            key = key.split("|")
            if key[6] == "GENE":
                k="|".join(key[0:6])
                exSize = int(key[7]); gSize = int(key[8])
                if key[0] in self.geneSeen:
                    EXP=self.geneExp[key[0]]
                    self.geneOut.write("%s %s %s | %s sense: %s %s %s antisense: %s %s %s\n" % (k,exSize,gSize,self.refType,EXP[0],EXP[1],EXP[2],EXP[3],EXP[4],EXP[5]))
                else:
                    self.geneOut.write("%s %s %s | %s sense: 0 0 0 antisense: 0 0 0\n" % (k,exSize,gSize,self.refType))
        







    def writeNovelGenes(self):

        if self.format == "SAM":
            self.spliceOut = open(self.prefix+'_splice.cnts','w')
            for s in self.spliceExp.keys():
                self.spliceOut.write("%s NOVEL %s\n" % (s,self.spliceExp[s]))
        else:
            self.minNovelCnt = 20 
            self.novelOut = open(self.prefix+'_novelAreas.cnts','w')
            self.passAreas = {}
            for chr in self.novelAreas:
                self.passAreas[chr]=dd(list)
                for area in self.novelAreas[chr]:
                    if self.novelAreas[chr][area] > self.minNovelCnt:
                        self.passAreas[chr][area] = []

            self.novelAreas.clear()

            self.mapHandle = open(self.mapFile)
            for line in self.mapHandle:
                line=line.split()
                myChr= line[2]; myPos = int(line[3]); myRound=myPos - myPos%1000; myStrand = line[5]; myMaps = int(line[7])
               
                if myMaps == 1 and myRound in self.passAreas[myChr]:
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
                    self.novelOut.write('%s %s %s +/- %s %s SPOTS: %s\n' % (chr,area,len(locs),pos,neg," ".join([str(l) for l in locs])))
                







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

        self.mapHandle.close()

        #self.fname.close()
        if self.samOut: self.samOut.close()
        if self.mapOut: self.mapOut.close()
        if self.AntiOut: self.AntiOut.close()
        if self.geneOut: self.geneOut.close()
        if self.spliceOut: self.spliceOut.close()
        if self.featOut: self.featOut.close()








