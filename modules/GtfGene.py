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
            self.seq = []
            self.seqIndex=0
            self.flankSeqs = None
            self.jPos    = 0 
            self.mutationKey = None
            
##########################################################################################################################################

#1)  ADDING ANOTHER LINE FROM A GTF FILE [ NOTICE YOU MUST START WITH GENE, THEN TRANSCRIPT OR EXON OK  


    def addGtfLine(self, line):
        if line.type=='transcript':
            self.transcripts.append([line.tranID,[]])

        elif line.type == 'exon':
            self.transcripts[-1][1].append((line.start,line.end))
            self.exons.append((line.start,line.end))
        else:
            return
  

###########################################################################################################################################

    def findMajorExons(self):
        tmpExons=[]
        if self.strand=="+":
            self.exons.sort()
            exStart=self.exons[0][0]; exEnd=self.exons[0][1]
            for i in range(1,len(self.exons)):
                if self.exons[i][0] <= exEnd:
                    if self.exons[i][1] > exEnd: exEnd = self.exons[i][1]
                else:
                    tmpExons.append((exStart,exEnd))
                    exStart=self.exons[i][0]; exEnd=self.exons[i][1]
            tmpExons.append((exStart,exEnd))
        else:
            self.exons.sort(reverse=True,key= lambda student: student [1])
            exStart=self.exons[0][0]; exEnd=self.exons[0][1]
            for i in range(1,len(self.exons)):
                if self.exons[i][1] >= exStart:
                    if self.exons[i][0] < exStart: exStart = self.exons[i][0]
                else:
                    tmpExons.append((exStart,exEnd))
                    exStart=self.exons[i][0]; exEnd=self.exons[i][1]
            tmpExons.append((exStart,exEnd))
        self.exons = tmpExons


#######################################################################################################################################################


    def findTranscriptOffsets(self):
        tmpTrans=[]
        if self.strand=="+":
            for t in self.transcripts:
                tOffsets=[]; gOffsets=[]; hOffsets=[]; k=0
                for loc in t[1]:
                    hOffsets.append([loc[0],loc[1]])
                    gOffsets.append([loc[0]-self.start,(loc[0]-self.start)+(loc[1]-loc[0])])
                    tOffsets.append([k,k+loc[1]-loc[0]])
                    k+=loc[1]-loc[0]+1
                tmpTrans.append([t[0],[tOffsets,gOffsets,hOffsets]])
        
        else:
            for t in self.transcripts:
                tOffsets=[]; gOffsets=[]; hOffsets=[]; k=0
                for loc in t[1]:
                    hOffsets.append([loc[1],loc[0]])
                    gOffsets.append([self.end-loc[1],self.end-loc[0]])
                    tOffsets.append([k,k+loc[1]-loc[0]])
                    k+=loc[1]-loc[0]+1
                tmpTrans.append([t[0],[tOffsets,gOffsets,hOffsets]])
        self.transcripts=tmpTrans


##############################################################################################################################################################

    def findFlankOffsets(self):

        self.flanks=[]
        if self.strand=="+":
            ## FRONT FLANK ## 
            tmpT=[[0,self.flankLen-1],[self.flankLen,self.flankLen+self.extend]]
            tmpH=[[self.start-self.flankLen,self.start-1],[self.start,self.start+self.extend]]
            tmpG=[["N"+str(self.flankLen),"N1"],[0,self.extend]]
            self.flanks.append(["5P",[tmpT,tmpG,tmpH]])
            
            ## BACK FLANK  ##
            tmpT=[[0,self.extend],[self.extend+1,self.flankLen+self.extend]]
            tmpH=[[self.end-self.extend,self.end],[self.end+1,self.end+self.flankLen]]
            tmpG=[[(self.end-self.start)-self.extend,(self.end-self.start)],[(self.end-self.start)+1,(self.end-self.start)+self.flankLen]]
            
            self.flanks.append(["3P",[tmpT,tmpG,tmpH]])
        else:
            ## FRONT FLANK ### 
            tmpT=[[0,self.flankLen-1],[self.flankLen,self.flankLen+self.extend]]
            tmpH=[[self.end+self.flankLen,self.end+1],[self.end,self.end-self.extend]]
            tmpG=[["N"+str(self.flankLen),"N1"],[0,self.extend]]
            self.flanks.append(["5P",[tmpT,tmpG,tmpH]])
            ## BACK FLANK  ##
            tmpT=[[0,self.extend],[self.extend+1,self.flankLen+self.extend]]
            tmpH=[[self.start+self.extend,self.start],[self.start-1,self.start-self.flankLen ]]
            tmpG=[[(self.end-self.start)-self.extend,(self.end-self.start)],[(self.end-self.start)+1,(self.end-self.start)+self.flankLen]]
            
            self.flanks.append(["3P",[tmpT,tmpG,tmpH]])
            


#######################################################################################################################################################################


    def findExonicOffsets(self):

        tmpT=[]; tmpG=[]; tmpH=[]; k=0; n=0;
        if self.strand=="+":
                
            for e in self.exons:
                tmpH.append([e[0],e[1]]); tmpG.append([e[0]-self.start,(e[0]-self.start)+(e[1]-e[0])]); tmpT.append([k,k+e[1]-e[0]])
                k+=e[1]-e[0]+1; n+=1
                if n<len(self.exons):
                    cut1=e[1]+1; cut2=self.exons[n][0]-1; span=cut2-cut1

                    if cut1-self.extend >= self.start and cut2+self.extend <= self.end:     EXT=self.extend
                    else:                                                                   EXT = min(cut1-self.start,self.end-cut2)

                    intT = [[0,EXT-1],[EXT,span+EXT],[span+EXT+1,span+EXT+EXT]]
                    intH = [[cut1-EXT,cut1-1],[cut1,cut2],[cut2+1,cut2+EXT]]
                    intG = [[cut1-EXT-self.start,cut1-1-self.start],[cut1-self.start,cut2-self.start],[cut2+1-self.start,cut2+EXT-self.start]]
                    self.introns.append([n,[intT,intG,intH]])
        else:
            for e in self.exons:
                tmpH.append([e[1],e[0]]);     tmpG.append([self.end-e[1],self.end-e[0]]);      tmpT.append([k,k+e[1]-e[0]])
                k+=e[1]-e[0]+1; n+=1
                if n<len(self.exons):
                    cut1=e[0]-1; cut2=self.exons[n][1]+1; span=cut1-cut2

                    if cut1+self.extend <= self.end and  cut2-self.extend >= self.start:    EXT=self.extend
                    else:                                                                   EXT=min(self.end-cut1, cut2 -self.start)
                    
                    intT = [[0,EXT-1],[EXT,span+EXT],[span+EXT+1,span+EXT+EXT]]
                    intH = [[cut1+EXT,cut1+1],[cut1,cut2],[cut2-1,cut2-EXT]]
                    intG = [[self.end-(cut1+EXT),self.end-(cut1+1)],[self.end-cut1,self.end-cut2],[self.end-(cut2-1),self.end-(cut2-EXT)]]
                    self.introns.append([n,[intT,intG,intH]])
        self.exons = [tmpT,tmpG,tmpH]
                    
        

    def validOffsets(self):

        if self.name==None:
            return False
        else:
            self.findMajorExons()
            self.findTranscriptOffsets()
            self.findFlankOffsets()
            self.findExonicOffsets()
            return True



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

        baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.seq = chrSeq[self.start-1:self.end]

        if self.strand == '-':
            self.seq.reverse()
            for i in range(len(self.seq)):
                self.seq[i] = baseComplement[self.seq[i].capitalize()]
        else:
            for i in range(len(self.seq)):
                self.seq[i] = self.seq[i].capitalize()

        ##########################################################################
        self.flankSeqs = [chrSeq[self.start-self.flankLen-1:self.start+self.extend],  chrSeq[self.end - self.extend-1:self.end+self.flankLen] ]

        if self.strand == '-':
            self.flankSeqs.reverse()
            for i in range(len(self.flankSeqs)):
                self.flankSeqs[i].reverse()
                for j in range(len(self.flankSeqs[i])):
                    self.flankSeqs[i][j] = baseComplement[self.flankSeqs[i][j].capitalize()]
        else:
            for i in range(len(self.flankSeqs)):
                for j in range(len(self.flankSeqs[i])):
                    self.flankSeqs[i][j] = self.flankSeqs[i][j].capitalize()

        
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




