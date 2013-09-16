#!/usr/bin/env python


import sys
from random import randrange
from random import choice
from collections import defaultdict as dd

##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class GtfGene:
    def __init__(self,line,readlen=100,flankLen=200,minOverlap=5):
        if line==None:
            self.name = None
        elif line=="NULL":
            self.start=0; self.end=0; self.name="NULL"
        elif line.type != 'gene':
            print line.type,"YO"
            print line.geneID, line.geneType,line.hugoName
            print "BAD INIT"
            sys.exit()
        else:

            self.name, self.chr, self.start, self.end, self.strand, self.type, self.status,self.hugo = line.geneID, line.chr, line.start,line.end,line.strand,line.geneType,line.geneStatus,line.hugoName 
            self.readlen = readlen; self.minOverlap = minOverlap; self.flankLen = flankLen;  self.extend = readlen-minOverlap; self.length = self.end-self.start+1
            
            if self.end - self.start <= self.extend:
                self.extend = (self.end - self.start)-1

            self.extendLen = self.extend
            self.transcripts, self.exons, self.introns, self.flanks = [],[],[],[]
            self.EXONS = set([])
            self.seq = []
            self.seqIndex=0
            self.flankSeqs = None
            self.jPos    = 0
            self.maxJump = 10
            self.mutationKey = None
            
###############################print###########################################################################################################

#1)  ADDING ANOTHER LINE FROM A GTF FILE [ NOTICE YOU MUST START WITH GENE, THEN TRANSCRIPT OR EXON OK  


    def addGtfLine(self, line):
        if line.type=='transcript':
            self.transcripts.append([line.tranID,[]])
        elif line.type == 'exon':
            
            
            if line.start < self.start:
                self.GENOME_RELATIVE = False

            self.transcripts[-1][1].append((line.start,line.end))
            self.EXONS.add((line.start,line.end))
        else:
            return
  

###########################################################################################################################################



    def checkOffsets(self):

        self.exons = sorted(list(self.EXONS))
        
        if self.exons[0][0] < self.start:
            for i in range(len(self.exons)):
                self.exons[i] = (self.exons[i][0] + self.start -1, self.exons[i][1] + self.start -1)
            for i in range(len(self.transcripts)):
                self.transcripts[i][1].sort()
                for j in range(len(self.transcripts[i][1])):
                    self.transcripts[i][1][j] = (self.transcripts[i][1][j][0] + self.start -1, self.transcripts[i][1][j][1] + self.start -1)
        if self.exons[-1][1] > self.end:
            self.end = self.exons[-1][1]




    def findMajorExons(self):
        majorExons = []; k=0; OPEN_FRAME = True
        exStart=self.exons[0][0]; exEnd=self.exons[0][1]
        for i in range(1,len(self.exons)):
            if self.exons[i][0] <= exEnd:
                if self.exons[i][1] > exEnd: exEnd = self.exons[i][1]
            else:
                majorExons.append((exStart,exEnd))
                exStart=self.exons[i][0]; exEnd=self.exons[i][1]
        
        majorExons.append((exStart,exEnd))
        self.exons = sorted(majorExons)
        majorSplice = []; myDists = []; tmpJxns = []

        while k < len(self.exons):
            if (self.exons[k][1]-self.exons[k][0]) > self.readlen:
                myDists.append((self.readlen))

                if majorSplice==[]:
                    majorSplice.append((self.exons[k][1]-(self.readlen),self.exons[k][1]))
                else:
                    majorSplice.append((self.exons[k][0],self.exons[k][0]+(self.readlen)))
            else:
                majorSplice.append(self.exons[k])
                myDists.append((self.exons[k][1]-self.exons[k][0]+1))
                    
            if sum(myDists[0:len(myDists)-1]) >= self.readlen and sum(myDists[1:len(myDists)]) >= self.readlen:
                tmpJxns.append(tuple(majorSplice))
                OPEN_FRAME = False
                while sum(myDists[1::]) >= self.readlen:
                    majorSplice.remove(majorSplice[0])
                    myDists.remove(myDists[0])
            k+=1

        if OPEN_FRAME and len(majorSplice)>1:
            tmpJxns.append(tuple(majorSplice))
        
        self.knownJunks  = sorted(tmpJxns)

            




#######################################################################################################################################################

    def findKnownJunctions(self):
        tmpJxns=set([]); myJumps=set([])
        for t in self.transcripts:
            k=0; OPEN_FRAME = True 
            tJxns = sorted(t[1]); mySplice=[]; myDists =[]
            while k < len(tJxns):

                if (tJxns[k][1]-tJxns[k][0]) > self.readlen:
                    myDists.append((self.readlen))
                    if mySplice==[]:
                        mySplice.append((tJxns[k][1]-(self.readlen),tJxns[k][1]))
                    else:
                        mySplice.append((tJxns[k][0],tJxns[k][0]+(self.readlen)))
                else:
                    mySplice.append(tJxns[k])
                    myDists.append((tJxns[k][1]-tJxns[k][0]+1))
                    
                if sum(myDists[0:len(myDists)-1]) >= self.readlen and sum(myDists[1:len(myDists)]) >= self.readlen:

                    if tuple(mySplice) not in self.knownJunks:
                        tmpJxns.add(tuple(mySplice))
                        for i in range(len(mySplice)-1):
                            myJumps.add((mySplice[i][1],mySplice[i+1][0]))
                    OPEN_FRAME = False
                    while sum(myDists[1::]) >= self.readlen:
                        mySplice.remove(mySplice[0])
                        myDists.remove(myDists[0])
                k+=1
            if OPEN_FRAME and len(mySplice)>1:
                if tuple(mySplice) not in self.knownJunks:
                    tmpJxns.add(tuple(mySplice))
                    for i in range(len(mySplice)-1):
                        myJumps.add((mySplice[i][1],mySplice[i+1][0]))
       
        self.knownJunks = sorted( self.knownJunks + [j for j in tmpJxns] )
        self.knownSplice = sorted([j for j in myJumps])



    def findNewJunctions(self):
        novelJxns=[]
        for i in range(len(self.exons)-2):
            for j in range(i+2,len(self.exons)):
                if j-i > self.maxJump: break
                sCand = (self.exons[i][1],self.exons[j][0])
                if sCand in self.knownSplice:
                    continue
                else:
                    if self.exons[i][1]-self.exons[i][0] > self.readlen:
                        firstEx = (self.exons[i][1] - self.readlen,self.exons[i][1])
                    else:
                        firstEx = self.exons[i]
                    if self.exons[j][1]-self.exons[j][0] > self.readlen:
                        secondEx = (self.exons[j][0],self.exons[j][0] + self.readlen)
                    else:
                        secondEx = self.exons[j]
                     
                    novelJxns.append((firstEx,secondEx))
        self.novelJxns = novelJxns
                
            
    def findIntrons(self):
        introns=[]
        for i in range(len(self.exons)):
            if i > 0:


                EX=self.exons[i]; PREV=self.exons[i-1]
                
                
                
                leftFlank  = max(PREV[1]-self.readlen,self.start)
                rightFlank = min(EX[0] + self.readlen, self.end)
                introns.append([(leftFlank,PREV[1]), (PREV[1]+1,EX[0]-1), (EX[0],rightFlank)])
               # introns.append([(PREV[1]-self.readlen,PREV[1]), (PREV[1]+1,EX[0]-1), (EX[0],EX[0]+self.readlen)])
        self.introns = introns



#######################################################################################################################################################################


                    

    def validOffsets(self):

        if self.name==None:
            return False
        else:
            self.checkOffsets()
            self.findMajorExons()
            self.findKnownJunctions()
            self.findNewJunctions()
            self.findIntrons() 
            return True






#######################################################################################################################################################################








































    def findSplicingInfo(self):
        self.spliceInfo=[]; self.jPos =0; self.codonOffset = 0
        for i in range(len(self.exons[1])):
            self.spliceInfo.append([self.exons[1][i],self.exons[2][i],set([(self.exons[1][i][0],self.exons[2][i][0])]),set([(self.exons[1][i][1],self.exons[2][i][1])])])

        tmpStarts=set([]); tmpEnds=set([])
        for t in self.transcripts:
            for i in range(len(t[1][1])):
                tmpStarts.add((t[1][1][i][0],t[1][2][i][0]))
                tmpEnds.add((t[1][1][i][1],t[1][2][i][1]))
       
        x=0; y=0; tmpStarts=sorted(list(tmpStarts)); tmpEnds=sorted(list(tmpEnds))
        for s in self.spliceInfo:
            while x<len(tmpStarts) and tmpStarts[x][0] < s[0][1]:
                s[2].add(tmpStarts[x]); x+=1
                if x==len(tmpStarts): break
            while y<len(tmpEnds) and tmpEnds[y][0] <= s[0][1]:
                s[3].add(tmpEnds[y]); y+=1
                if y==len(tmpEnds): break



    def findJxns(self,gPos):
       
        if gPos < 0:
            if self.strand == "+": return "5P-FLNK","NA",['+',self.length,self.start]
            else:                  return "5P-FLNK","NA",['-',self.length,self.end]
        elif gPos > self.length:
            if self.strand == "+": return "3P-FLNK","NA",['+',self.length,self.end]
            else:                  return "3P-FLNK","NA",['-',self.length,self.start]

        
        while self.jPos < len(self.spliceInfo):
            
            
              
            if self.spliceInfo[self.jPos][0][0] <= gPos and gPos <= self.spliceInfo[self.jPos][0][1]:

                return "EXONIC", (self.codonOffset + (gPos - self.spliceInfo[self.jPos][0][0])) %3,  self.spliceInfo[self.jPos]

            elif gPos < self.spliceInfo[self.jPos][0][0]:

                if self.jPos == 0:
                    jHere = self.spliceInfo[self.jPos]
                    return "INTRONIC",'NA',[[jHere[0][0],jHere[0][0]],[jHere[1][0],jHere[1][0]],"NULL"]
                elif self.spliceInfo[self.jPos-1][0][1] < gPos:

                    jPrev = self.spliceInfo[self.jPos-1]; jNext = self.spliceInfo[self.jPos]
                
                    return "INTRONIC", 'NA',[[jPrev[0][1],jNext[0][0]],[jPrev[1][1],jNext[1][0]],"NULL"]
             


            else:
                self.jPos+=1
                self.codonOffset += ( self.spliceInfo[self.jPos][0][1] - self.spliceInfo[self.jPos][0][0] + 1 ) % 3 




    def getSeqFromChr(self,chrSeq):

        self.seq = [base.capitalize() for base in chrSeq[self.start-1:self.end]]

        ##########################################################################
        
        
###########################################################################################################################################
    
    def mutateSeq(self,exRate,intRate):

        if not self.seq:
            print "NO SEQ YET"
            sys.exit()

        mutList = dd(list); baseList = ["A","C","G","T"]

        mutation_rates = [0.15, 0.25, 0.35, 0.45, 0.5, 0.55, 0.65, 0.75, 0.9, 1.0, 1.0, 1.0 , 1.0, 1.0 ]
        for i in range(len(self.exons[1])):
            start=self.exons[1][i][0]; end = self.exons[1][i][1]; span=end-start; k= 0 
            
            if span < 10: continue
            
            while (k+1.0) / span < exRate:

                gPos = randrange(start,end)
                myRef = self.seq[gPos]
                if mutList[gPos] or myRef == "N":
                    continue
                else:
                    if self.strand == "+":  hPos = self.exons[2][i][0] + ( gPos - start )
                    else:                   hPos = self.exons[2][i][0] - ( gPos - start )
                    
                    baseList.remove(myRef)
                    mutList[gPos] = [ "EXONIC",gPos, hPos, myRef, choice(baseList), choice(mutation_rates),0,0 ]
                    self.seq[gPos] = 'x'
                    baseList.append(myRef)
                    k+=1


        for i in range(len(self.introns)):
            start = self.introns[i][1][1][1][0]; end = self.introns[i][1][1][1][1]; span = end-start; k=0
            
            if span < 10: continue
            
            while (k+1.0) / span < intRate:
                
                gPos = randrange(start,end)
                myRef = self.seq[gPos]

                if mutList[gPos] or myRef == "N":
                    continue
                else:
                    if self.strand == "+":  hPos = self.introns[i][1][2][1][0] + ( gPos - start)
                    else:                   hPos = self.introns[i][1][2][1][0] - ( gPos - start)
                    
                    baseList.remove(myRef)
                    mutList[gPos] = [ "INTRONIC",gPos, hPos, myRef, choice(baseList), choice(mutation_rates),0,0 ]
                    self.seq[gPos] = 'x'
                    baseList.append(myRef)
                    k+=1
        self.mutationKey =  mutList




