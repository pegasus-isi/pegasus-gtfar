#!/usr/bin/env python


import sys
#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/PIPELINE_v2/gtFar/python_src/gtfar_modules')
#from tools.gtTools import *



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
            print "BAD INIT"
            sys.exit()
        else:
            self.name, self.chr, self.start, self.end, self.strand, self.type, self.status,self.hugo = line.geneID, line.chr, line.start,line.end,line.strand,line.geneType,line.geneStatus,line.hugoName
            
            self.readlen = readlen; self.minOverlap = minOverlap; self.flankLen = flankLen;  self.extend = readlen-minOverlap; self.length = self.end-self.start+1
            
            if self.end - self.start <= self.extend:
                self.extend = (self.end - self.start)-1


            self.transcripts, self.exons, self.introns, self.flanks = [],[],[],[]
            self.seq = []
            self.seqIndex=0
            self.flankSeqs = None
            
            
##########################################################################################################################################

#1)  ADDING ANOTHER LINE FROM A GTF FILE 


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
                if self.exons[i][0] < exEnd:
                    if self.exons[i][1] > exEnd: exEnd = self.exons[i][1]
                else:
                    tmpExons.append((exStart,exEnd))
                    exStart=self.exons[i][0]; exEnd=self.exons[i][1]
            tmpExons.append((exStart,exEnd))
        else:
            self.exons.sort(reverse=True,key= lambda student: student [1])
            exStart=self.exons[0][0]; exEnd=self.exons[0][1]
            for i in range(1,len(self.exons)):
                if self.exons[i][1] > exStart:
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
        self.spliceInfo=[]
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
            





    def getSeqFromChr(self,chrSeq,index):

        baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.seq = chrSeq[self.start-1:self.end]
        self.seqIndex = index

        if self.strand == '-':
            self.seq.reverse()
            for i in range(len(self.seq)):
                self.seq[i] = baseComplement[self.seq[i].capitalize()]
        else:
            for i in range(len(self.seq)):
                self.seq[i] = self.seq[i].capitalize()
        if index == 1:
            self.seq = ['N']+self.seq

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

        

    def forceIndex(self,idx):
        if self.seqIndex != idx:
            if self.seqIndex == 1 and idx == 0:
                self.seq = self.seq[1::]
                self.seqIndex = 0
            elif self.seqIndex == 0 and idx == 1:
                self.seq = ['N']+self.seq
                self.seqIndex = 1







