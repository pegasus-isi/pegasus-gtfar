#!/usr/bin/env python

import os
import sys
import difflib
from collections import defaultdict as dd
from math import fabs

from Sequence import *

#from ..gtTools.seq import *


##########################################################################################################################################
#####################################################  MAPPED READ  CLASS START ############################################################
##########################################################################################################################################



class MapRead:
    def __init__(self,type,sense,ID,qual,subs,tMaps,readLines):

        self.ID = ID; self.sense = sense; self.qual = qual; self.subs = subs; self.features = tMaps; self.readLines = readLines; self.PALINDROMIC = False; self.TRIMMED = False;

        self.samStrings = []; self.locStrings = []
       

        if type == 'MAPPING' or type == 'SAM':
        
            if len(self.readLines) > 1:  self.disambiguate()
            
            self.hgLocs = [self.readLines[0][0]]; self.geneLocs = [[self.readLines[0][1]]]; self.seqs = [self.readLines[0][2]]
            
            for i in range(1,len(self.readLines)):

                if self.readLines[i][0] != self.hgLocs[-1]:
                    self.hgLocs.append(self.readLines[i][0])
                    self.seqs.append(self.readLines[i][2])
                    self.geneLocs.append([self.readLines[i][1]])
                else:
                    if self.readLines[i][1] != self.geneLocs[-1][-1]:
                        self.geneLocs[-1].append(self.readLines[i][1])
             
            if len(self.hgLocs) == 1:
                self.uniq = True; self.multi=False
                if len(self.geneLocs[0])!=1:
                    self.multi=True
            else:
                self.uniq=False; self.multi=True

###################################################################################################################################



                    


#####################################################################################################################################################

    def palindromic(self):

            ### CHECK FOR PALINDROMIC AMBIGUITY ###

        if self.readLines[0][0] != self.readLines[-1][0] and self.readLines[0][0][0] == self.readLines[-1][0][0]:

            if self.readLines[0][0][2][-1::-1] == self.readLines[-1][0][2]:

                n=0; self.PALINDROMIC = True
                while n < len(self.readLines):
                    if self.readLines[n][0][1] != "+":
                        self.readLines.remove(self.readLines[n])
                    else:
                        n+=1
            else:
                return


    def startTrim(self,n,minLen,maxLen):
        seqLen = 0; myStrand = self.readLines[0][n][1]; newSpots = []; readLen = len(self.readLines[0][2][0])
        for i in range(0,minLen,2):
            x1 = set([]); x2=set([])
            for j in range(len(self.readLines)):
                x1.add(self.readLines[j][n][2][i]);  x2.add(self.readLines[j][n][2][i+1])
            if len(x1) == 1:
                x1 = x1.pop()
                if len(x2) == 1:    x2=x2.pop()
                else:
                    if x1 < min(x2): x2 = min(x2)
                    else:            x2 = max(x2)
                newSpots.extend([x1,x2])
                seqLen += fabs(x2-x1)+1
        if myStrand == "+":
            return newSpots,(0,seqLen)
        else:
            return newSpots,(readLen-seqLen,readLen)
                
                
    def tailTrim(self,n,minLen,maxLen):
        seqLen = 0; myStrand = self.readLines[0][n][1]; newSpots = []; readLen = len(self.readLines[0][2][0])
        for i in range(0,minLen,2):
            x1 = set([]); x2=set([])
            for j in range(len(self.readLines)):
                sData=self.readLines[j][n][2]
                x1.add(sData[len(sData)-(i+2)]); x2.add(sData[len(sData)-(i+1)])

            if len(x2) == 1:
                x2 = x2.pop()
                if len(x1) == 1:    x1=x1.pop()
                else:
                    if x2 > max(x1):     x1=max(x1)
                    else:           x1=min(x1)
                newSpots.extend([x1,x2])
                seqLen += fabs(x2-x1)+1

        if myStrand == "+":
            return newSpots,(readLen-seqLen,readLen)
        else:
            return newSpots,(0,seqLen)

        



    def trim(self):

        hgEnds   = set([h[0][2][-1] for h in self.readLines]);
        hgStarts = set([h[0][2][0]  for h in self.readLines]);
        minLen   = min([len(h[0][2]) for h in self.readLines])
        maxLen   = max([len(h[0][2]) for h in self.readLines])
       
        geEnds   = set([h[1][2][-1] for h in self.readLines]);
        geStarts = set([h[1][2][0]  for h in self.readLines]);
       

        if len(hgStarts) > 1 and len(hgEnds) > 1:
            return 
        else:
            if len(hgStarts) == 1:   hgResult = self.startTrim(0,minLen,maxLen)
            elif len(hgEnds) == 1:   hgResult = self.tailTrim(0,minLen,maxLen)
            else:                    return 
            if len(geStarts) == 1:   geResult = self.startTrim(1,minLen,maxLen)
            elif len(geEnds) == 1:   geResult = self.tailTrim(1,minLen,maxLen)
            else:                    return
        
        if len(hgResult[0]) == len(geResult[0]) and hgResult[1] == geResult[1]:
            self.TRIMMED=True
            fd = self.readLines[0]; A=int(hgResult[1][0]); B=int(hgResult[1][1])
            self.readLines = [[(fd[0][0],fd[0][1],hgResult[0]), (fd[1][0],fd[1][1],geResult[0]), (fd[2][0][A:B],fd[2][1][A:B])]]










    def trimLeft(self,hgLens,maxMap,minMap):
        if maxMap[0] == minMap[0]:

            newChr = []; newGene = []; self.TRIMMED=True
            for i in range(len(self.readLines[0][0][2])):
                tmpHg=set([]); tmpGene=set([])

                for j in range(len(self.readLines)):

                    tmpHg.add(self.readLines[j][0][2][i])
                    tmpGene.add(self.readLines[j][1][2][i])
                    
                if len(tmpHg) == 1:
                    newChr.append(tmpHg.pop())
                    newGene.append(tmpGene.pop())
                for j in range(len(self.readLines)):

                    tmpGene.add(self.readLines[j][1][2][i])
                if len(tmpGene) == 1: newGene.append(tmpGene.pop())

            if len(newGene) == 1:
                addedCands = sorted([ (self.readLines[j][1][2][1],self.readLines[j][0][2][1]) for j in range(len(self.readLines)) ])
                newGene.append(addedCands[0][0])
                newChr.append(addedCands[0][1])
            seqDist=0;
            for i in range(0,len(newGene),2):
                seqDist+=(newGene[i+1]-newGene[i])+1
            self.qual = self.qual[0:seqDist]
            for i in range(len(self.readLines)):
                self.readLines[i][0] = (self.readLines[i][0][0],self.readLines[i][0][1],newChr)
                self.readLines[i][1] = (self.readLines[i][1][0],self.readLines[i][1][1],newGene)
                self.readLines[i][2] = (self.readLines[i][2][0][0:seqDist], self.readLines[i][2][1][0:seqDist])
        else:
        
            ### TEST IT OUT ###

            newMap  = [self.readLines[0][0][0],self.readLines[0][0][1],maxMap[1][0:minMap[0]]];
            newGene = [self.readLines[0][1][0],self.readLines[0][1][1],maxMap[2][0:minMap[0]]];
            NEWLEN=int(sum([fabs(newMap[2][i+1]-newMap[2][i])+1 for i in range(0,len(newMap[2]),2)]))
            self.qual = self.qual[0:NEWLEN]
            
            for i in range(len(self.readLines)):
                self.readLines[i][0] = newMap
                self.readLines[i][1] = newGene
                self.readLines[i][2] = ( self.readLines[i][2][0][0:NEWLEN], self.readLines[i][2][1][0:NEWLEN] )


    def trimRight(self,hgLens,maxMap,minMap):
        if maxMap[0] == minMap[0]:
            newChr = []; newGene = []; self.TRIMMED=True; k=0
            for i in range(len(self.readLines[0][0][2])):
                tmpHg=set([]); tmpGene=set([])
                for j in range(len(self.readLines)):
                    tmpHg.add(self.readLines[j][0][2][i])
                    tmpGene.add(self.readLines[j][1][2][i])
                    
                if len(tmpHg) == 1:
                    newChr.append(tmpHg.pop())
                    newGene.append(tmpGene.pop())
                else:
                    k=i


            if len(newGene) == 1:
                addedCands = sorted([ (self.readLines[j][1][2][i],self.readLines[j][0][2][i]) for j in range(len(self.readLines)) ]) 
                newGene = [addedCands[0][0]] + newGene
                newChr  = [addedCands[0][1]] + newChr


            seqDist=0
            for i in range(0,len(newGene),2):
                seqDist+=(newGene[i+1]-newGene[i])+1

            seqStart = len(self.qual)-seqDist
            self.qual = self.qual[seqStart::]
            for i in range(len(self.readLines)):
                self.readLines[i][0] = (self.readLines[i][0][0],self.readLines[i][0][1],newChr)
                self.readLines[i][1] = (self.readLines[i][1][0],self.readLines[i][1][1],newGene)
                self.readLines[i][2] = (self.readLines[i][2][0][seqStart::], self.readLines[i][2][1][seqStart::])
        

        else:
            TRIM = maxMap[1][0:maxMap[0]-minMap[0]];    self.TRIMMED=True; TRIMLEN = int(sum([fabs(TRIM[i+1]-TRIM[i])+1 for i in range(0,len(TRIM),2)])) 
            newMap  = [self.readLines[0][0][0],self.readLines[0][0][1],maxMap[1][(maxMap[0]-minMap[0])::]];  self.qual = self.qual[TRIMLEN::]
               
            for i in range(len(self.readLines)):

                self.readLines[i][0] = newMap; gmap = self.readLines[i][1]
                        
                if len(self.readLines[i][1][2]) > minMap[0]:    self.readLines[i][1] = (gmap[0],gmap[1],gmap[2][len(gmap[2])-len(newMap[2])::])
                else:                                           self.readLines[i][1] = (gmap[0],gmap[1],[gmap[2][0]+TRIMLEN]+gmap[2][1::])
                        
                self.readLines[i][2] = ( self.readLines[i][2][0][TRIMLEN::], self.readLines[i][2][1][TRIMLEN::] )
                            



    def removeDuplicates(self):
        n=1; self.readLines.sort()
        while n < len(self.readLines):
            if self.readLines[n] == self.readLines[n-1]:
                self.readLines.remove(self.readLines[n])
            else:
                n+=1







    def disambiguate(self):

            ### CHECK FOR PALINDROMIC AMBIGUITY ###
        
        self.readLines.sort()


        chrSet=set([]); geneSet=set([]); refSet=set([]); cStrand=set([])
        for i in range(len(self.readLines)):
            chrSet.add(self.readLines[i][0][0]); geneSet.add(self.readLines[i][1][0]); refSet.add(self.readLines[i][2][1]); cStrand.add(self.readLines[i][0][1]) #;  gStrand=self.readLines[i][1][1]; 

        if len(chrSet) == 1:


            if self.readLines[0][0][2][-1::-1] == self.readLines[-1][0][2]:
                self.palindromic()


            elif len(geneSet) == 1 and len(refSet)==1 and len(cStrand) == 1:

                #READ=self.readLines[0][2][0]; REF=refSet.pop()
                #cStrand = self.readLines[0][0][1]; rStrand=self.readLines[0][1][1]
                self.trim()
                hgEnds   = set([h[0][2][-1] for h in self.readLines]);
                hgStarts = set([h[0][2][0]  for h in self.readLines]);
                hgLens   = [(len(h[0][2]),h[0][2],h[1][2]) for h in self.readLines];
                maxMap   = max(hgLens);
                minMap   = min(hgLens);

            #    if len(hgStarts) == 1:
             #       self.trimLeft(hgLens,maxMap,minMap)

              #  elif len(hgEnds) == 1:
               #     self.trimRight(hgLens,maxMap,minMap)
            
            self.removeDuplicates()

                        
                
















##################################################################################################################################################################################
##################################################################################################################################################################################
#########################################################################  READ PROCESSES ########################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################






##################################################################################################################################################################################
#########################################################################  READ PROCESSES ########################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################





    def visualizeStrings(self):

        chrDict ={'chr1': 1, 'chr2': 2, 'chr3': 3,'chr4': 4, 'chr5': 5, 'chr6': 6,'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr10': 10, 'chr11': 11, 'chr12': 12, 
                'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18,'chr19': 19, 'chr20': 20, 'chr21': 21,'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrMT': 25}
 
        if  self.uniq:
            samScr='255'; hgType="UNIQ"
        else:
            samScr='1';   hgType="AMBIG"

        cigarStr = ''; spotStr = ''

        for i in range(len(self.hgLocs)):
            hg = self.hgLocs[i]; mySeq = self.seqs[i][0]
            if hg[1]=='+':
                samSeq = mySeq;  spots = hg[2]; samStrand = '0'
            else:
                samSeq = revComp(mySeq); spots = hg[2][-1::-1]; samStrand='16'; self.qual = self.qual[-1::-1]


            for k in range(0,len(spots),2):
                cigarStr+=str(spots[k+1]-spots[k]+1)+"M"
                if k+2< len(spots):
                    cigarStr+= str(spots[k+2] - spots[k+1] -1)+"N"

 #           print self.hgLocs, self.ID,self.PALINDROMIC, self.TRIMMED 


            self.samStrings.append(self.ID+'\t'+samStrand+'\t'+hg[0]+'\t'+str(spots[0])+'\t'+samScr+'\t'+cigarStr+'\t*\t0\t0\t'+samSeq+'\t'+'\t'+self.qual+'\tNM:i:'+str(self.subs))
        
            for j in range(len(self.geneLocs[i])):
                if len(self.geneLocs[i]) > 1: gType = "MULTI"
                else:                         gType = "SINGLE"
                gene = self.geneLocs[i][j]
                if gene[1] == '+':
                    mySeq = self.seqs[i][0]; myRef = self.seqs[i][1]
                else:
                    mySeq = revComp(self.seqs[i][0]); myRef = revComp(self.seqs[i][1]); self.qual[-1::-1]

                spotStr=','.join([str(s) for s in gene[2]])
                geneString=" ".join([self.ID,hgType,gType,hg[0],hg[1],gene[0],gene[1],spotStr,"SEQ/REF/QUAL",mySeq,myRef,self.qual,str(self.subs),str(chrDict[hg[0]])])
                if geneString not in self.locStrings:   self.locStrings.append(geneString)
                
                 
#########################################################################################################################################################################################################


































































































                

                    



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






