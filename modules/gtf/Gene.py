#!/usr/bin/env python


import sys




##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class Gene:
    def __init__(self,line,readlen,flankLen=200,minOverlap=5):
        if line==None:
            self.name = None
        elif line.type != 'gene':
            print "BAD INIT"
            sys.exit()
        else:
            self.name, self.chr, self.start, self.end, self.strand, self.type, self.status,self.hugo = line.geneID, line.chr, line.start,line.end,line.strand,line.geneType,line.geneStatus,line.hugoName
            
            self.readlen = readlen; self.minOverlap = minOverlap; self.flankLen = flankLen;  self.extend = readlen-minOverlap
            
            self.transcripts, self.exons, self.introns, self.flanks = [],[],[],[]
            
            

    def add(self, line):
        if line.type=='transcript':
            self.transcripts.append([line.tranID,[]])

        elif line.type == 'exon':
            self.transcripts[-1][1].append((line.start,line.end))
            self.exons.append((line.start,line.end))
        else:
            return
    
    def finish(self):
        if self.name==None:
            return False
        
        tmpTrans=[]
        tmpExons=[]

        ## WE MAKE CATTED EXONS, INTRONIC SEQS, AND TRANSCRIPTS ## 

        ## 0) FIND MAJOR EXONS ## 

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

        ## 1) MAKING THE TRANSCRIPT SEQUENCES ## 
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

        ## 2) MAKING THE EXONIC CATS AND INTRONS ##
        tOffsets=[]; gOffsets=[]; hOffsets=[]; k=0;n=0
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
           
            for e in tmpExons:
                hOffsets.append([e[0],e[1]])
                gOffsets.append([e[0]-self.start,(e[0]-self.start)+(e[1]-e[0])])
                tOffsets.append([k,k+e[1]-e[0]])
                k+=e[1]-e[0]+1; n+=1
                if n<len(tmpExons):
                    cut1=e[1]+1; cut2=tmpExons[n][0]-1; span=cut2-cut1

                    if cut1-self.extend >= self.start and cut2+self.extend <= self.end:
                        EXT=self.extend
                    else:
                        EXT = min(cut1-self.start,self.end-cut2)


                    tmpT = [[0,EXT-1],[EXT,span+EXT],[span+EXT+1,span+EXT+EXT]]
                    tmpH = [[cut1-EXT,cut1-1],[cut1,cut2],[cut2+1,cut2+EXT]]
                    tmpG = [[cut1-EXT-self.start,cut1-1-self.start],[cut1-self.start,cut2-self.start],[cut2+1-self.start,cut2+EXT-self.start]]
                    self.introns.append([n,[tmpT,tmpG,tmpH]])
            self.exons = [tOffsets,gOffsets,hOffsets]
            
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
            

            for e in tmpExons:
                hOffsets.append([e[1],e[0]])
                gOffsets.append([self.end-e[1],self.end-e[0]])
                tOffsets.append([k,k+e[1]-e[0]])
                k+=e[1]-e[0]+1; n+=1
                if n<len(tmpExons):
                    cut1=e[0]-1; cut2=tmpExons[n][1]+1; span=cut1-cut2

                    if cut1+self.extend <= self.end and  cut2-self.extend >= self.start:
                        EXT=self.extend
                    else:
                        EXT=min(self.end-cut1, cut2 -self.start)
                    
                    tmpT = [[0,EXT-1],[EXT,span+EXT],[span+EXT+1,span+EXT+EXT]]
                    tmpH = [[cut1+EXT,cut1+1],[cut1,cut2],[cut2-1,cut2-EXT]]
                    tmpG = [[self.end-(cut1+EXT),self.end-(cut1+1)],[self.end-cut1,self.end-cut2],[self.end-(cut2-1),self.end-(cut2-EXT)]]
                    self.introns.append([n,[tmpT,tmpG,tmpH]])
            self.exons = [tOffsets,gOffsets,hOffsets]
                    
        return True

##########################################################################################################################################
#####################################################  GENE CLASS END  ###################################################################
##########################################################################################################################################


