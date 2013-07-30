#!/usr/bin/env python


import sys
import os 
import difflib

from Sequence import *
from Utilities import *
from GtfLine import *
from GtfGene import *
from random import randrange
from random import random 
from random import choice

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix,readlen,TYPE="ALL",printKEY=True):
        
        self.fName = open(fileHandle);  self.prefix = prefix;  self.readlen = readlen; self.PRINTKEY = printKEY
        
        self.index =0; self.genes = []; self.geneKey={}; self.seq = []
        
        self.open = True; self.FIRSTPRINT=True; self.notes = None; self.mutFile = None; self.simReads = None; self.counter = 0
        
        ##############################################


        tmpLine=self.fName.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr

############################################################################################################################################


    def loadGenesOnChromosome(self):
        while self.chr == self.line.chr:
            self.gene = GtfGene(self.line,self.readlen)
            while self.line.valid and self.line.geneID ==  self.gene.name:
                self.gene.addGtfLine(self.line)
                self.line = GtfLine(self.fName.readline())
            if self.gene.validOffsets():
                self.genes.append(self.gene)
                self.geneKey[self.gene.name] = self.gene 

             
    
    def startNextChromosome(self):

        if self.line.chr == 'NA':
            self.chr =  'NA'
            self.open = False
        else:
            self.chr = self.line.chr; self.genes=[]; self.seq=[]; self.geneKey={}


    def addFasta(self,filePath):
        c=open(filePath)
        if self.chr != c.readline().strip().split(">")[1]:
            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])

#####################################################################################################################################################3


    def uniquifySeqs(self,SILENT=False):
        if self.notes == None and SILENT != True:
            self.notes=open(self.prefix+'.notes','w')
        
        i=0; lastPos=GtfGene("NULL"); lastNeg=GtfGene("NULL"); removeList=[]; warnList=[]
        while i < len(self.genes):
            myGene=self.genes[i]
            if myGene.strand == "+":
                if myGene.start>lastPos.end:
                    lastPos=myGene; i+=1
                else:
                    if myGene.end<= lastPos.end:
                        removeList.append([myGene,lastPos,self.chr,'+']);  self.genes.remove(self.genes[i])
                    else:
                        warnList.append([myGene,lastPos,self.chr,'+']); i+=1; lastPos=myGene

            else:
                if myGene.start>lastNeg.end:
                    lastNeg=myGene; i+=1
                else:
                    if myGene.end<= lastNeg.end:
                        removeList.append([myGene,lastNeg,self.chr,'-']);  self.genes.remove(self.genes[i])
                    else:
                        warnList.append([myGene,lastNeg,self.chr,'-']); i+=1; lastPos=myGene
        if self.notes != None:
            for n in removeList:
                self.notes.write('Note: %s is a subsequence of %s and will not be annotated (%s %s [%s to %s] vs. [%s to %s])\n' % (n[0].name,n[1].name,n[2],n[3],n[0].start,n[0].end,n[1].start,n[1].end))
            for n in warnList:
                self.notes.write('Note: %s shares sequence with %s and may cause ambiguity (%s %s [%s to %s] vs. [%s to %s])\n' % (n[0].name,n[1].name,n[2],n[3],n[0].start,n[0].end,n[1].start,n[1].end))

            

###########################################################################################################################################################################################################



    def concatSeq(self,gSeq,restarts):
        mySeq=''
        for restart in restarts: mySeq+=gSeq[restart[0]:restart[1]+1]
        return mySeq


    def junctionLocations(self,jumps):
        
        locStrs=[[] for x in range(len(jumps))]
        for i in range(len(jumps)):
            for j in range(len(jumps[i])):
                locStrs[i].append(str(jumps[i][j][0])+"-"+str(jumps[i][j][1]))
            locStrs[i]=",".join(locStrs[i])
        return " ".join(locStrs)


########################################################################################################################################################
##################################### METHODS TO EVALUATE MUTATIONS  ###################################################################################
########################################################################################################################################################

# INPUT IS A DICTIONARY OF MUTATION KEYS - THE KEY IS THE GENE ID - THE VALUE IS LIST A MUTATION CANDIDATE OBJECTS #

    def evaluateMutations(self,mutations,minCov=3,minMut=0.5,minUrn=0.1):
        if not self.mutFile:    self.mutFile=open(self.prefix+'.mutations','w')
        for gene in self.genes:
            if len(mutations[gene.name]) > 0:
                gene.getSeqFromChr(self.seq)
                gene.findSplicingInfo()
                for cand in mutations[gene.name]:
                    ## OK GET THESE DISTS ##
                    cand.seqType,cand.codonOffset,candJxns = gene.findJxns(cand.pos)
                    cand.evalJxns(candJxns)
                    
                    if cand.seqType == "EXONIC":        cand.addRef(gene.seq[cand.pos-2:cand.pos+3])
                    elif cand.seqType == "INTRONIC":    cand.addRef(gene.seq[cand.pos])
                    else:
                        if cand.seqType ==  "5P-FLNK":
                                cand.addRef(gene.flankSeqs[0][gene.flankLen - cand.spliceDist])
                        elif cand.seqType == "3P-FLNK":
                                cand.addRef(gene.flankSeqs[1][gene.extendLen + cand.spliceDist  ])

                    cand.evalStats()

                    self.mutFile.write('%s:%s %s | %s %s:%s %s | ' % (gene.chr,cand.hgPos,gene.strand,gene.hugo,gene.name,cand.pos,cand.seqType))
                    self.mutFile.write('%s %s %s %g %s | ' % (cand.refBase,cand.altBase,cand.cov,cand.altRate,",".join([str(s) for s in cand.cnts])))
             

                    self.mutFile.write('CODON: %s %s | ' % ("".join(cand.codonSeq), cand.codonOffset))
                    if cand.spliceType=="BEG":  self.mutFile.write('j-dist: %s BEG %s %s %s | ' % (cand.spliceDist,cand.spliceClass,cand.start_tuple[0],cand.start_tuple[1]))
                    else:                       self.mutFile.write('j-dist: %s END %s %s %s | ' % (cand.spliceDist,cand.spliceClass,cand.end_tuple[0],cand.end_tuple[1]))
                    self.mutFile.write('uVal/qVal: %g %g %g %g | ' % (cand.urnRef,cand.urnAlt,cand.qualRef,cand.qualAlt))

                    if cand.cov < minCov:       self.mutFile.write('EDIT=FALSE (Low Cov) \n')
                    elif cand.altRate < minMut: self.mutFile.write('EDIT=FALSE (Low Mutation Rate) \n')
                    elif cand.urnAlt < minUrn:  self.mutFile.write('EDIT=FALSE (Redundnant Positions) \n')
                    else:
                        self.mutFile.write('EDIT=TRUE \n')
                        if cand.seqType != "5P-FLNK" and cand.seqType !="3P-FLNK" and cand.seqType != "FLANK":
                            gene.seq[cand.pos] = cand.altBase
                        else:
                            if candJxns[0] == '3P': gene.flankSeqs[1][gene.extendLen + cand.spliceDist] = cand.altBase 
                            

#######################################################################################################################################################
################################## METHODS TO PRINT OUT ANNOTATION AND SEQUENCES ######################################################################
#######################################################################################################################################################





    def printAnnotation(self,TYPE='ALL'):
        self.printTYPE=TYPE

        if self.FIRSTPRINT:
            self.FIRSTPRINT = False
            if TYPE=='ALL' or TYPE=='EXONIC':
                self.genePrint   = open(self.prefix+'_geneSeqs.fa','w')
                self.exonPrint   = open(self.prefix+'_exonSeqs.fa','w')
                self.catsOnly    = open(self.prefix+'_catsOnly.fa','w')
            if TYPE=='ALL' or TYPE=='INTRONIC':
                #self.genePrint   = open(self.prefix+'_geneSeqs.fa','w')
                self.intronPrint = open(self.prefix+'_intronSeqs.fa','w')
            
            if self.PRINTKEY:
                self.keyFile     = open(self.prefix+'.key','w')

        for gene in self.genes:
            if not gene.seq:
                gene.getSeqFromChr(self.seq)
            self.printSequence(gene)


    def printSequence(self,gene):

        geneSeq="".join(gene.seq)
        geneStr=gene.name+"|"+gene.hugo+"|"+gene.strand+"|"+gene.chr+"|"+gene.type+"|"+gene.status+"|"
        exCatSeq=self.concatSeq(geneSeq,gene.exons[1]); tCnt=0
        catString=geneStr+"GENE|"+str(len(exCatSeq))+"|"+str(len(geneSeq))+"|"

        if self.printTYPE != 'EXONIC':

            flankSeqs = ["".join(seq) for seq in gene.flankSeqs]; iCnt=0
            for i in range(len(gene.flankSeqs)):
                flankSeq="".join(gene.flankSeqs[i])
                flankStr=geneStr+"FLNK"+str(i)+":"+str(gene.flankLen)+"&"+str(gene.extend)+"|"
                if i==0:
                    self.intronPrint.write(">%s\n%s\n" % (flankStr+"5p",flankSeq))
                else:
                    self.intronPrint.write(">%s\n%s\n" % (flankStr+"3p",flankSeq))
                if self.PRINTKEY: self.keyFile.write("%s %s\n" % (flankStr,self.junctionLocations(gene.flanks[i][1])))
            
            for i in gene.introns:
                intSeq= self.concatSeq(geneSeq,i[1][1]); iCnt+=1
                intStr=geneStr+"ITRN"+str(iCnt)+":"+str(len(intSeq))+"|"
                self.intronPrint.write(">%s\n%s\n" % (intStr+"intRon",intSeq))
                if self.PRINTKEY: self.keyFile.write("%s %s\n" % (intStr,self.junctionLocations(i[1])))

        if self.printTYPE != 'INTRONIC':

            self.exonPrint.write(">%s\n%s\n" % (catString+"CATSEQ",exCatSeq))
            self.catsOnly.write(">%s\n%s\n" % (catString+"CATSEQ",exCatSeq))
            self.genePrint.write(">%s\n%s\n" % (catString+"FULLSEQ",geneSeq))
            if self.PRINTKEY:
                self.keyFile.write("%s %s\n" % (catString,self.junctionLocations(gene.exons)))
            for t in gene.transcripts:
                tranSeq= self.concatSeq(geneSeq,t[1][1]); tCnt+=1
                tranStr=geneStr+"TRAN"+str(tCnt)+":"+str(len(tranSeq))+"|"
                self.exonPrint.write(">%s\n%s\n" % (tranStr+t[0],tranSeq))
                if self.PRINTKEY: self.keyFile.write("%s %s\n" % (tranStr,self.junctionLocations(t[1])))
        
        

##################################################################################################################################################################
       
       


#######################################################################################################################################################
################################## METHODS TO SIMULATE READS FROM GTF FILE   ##########################################################################
#######################################################################################################################################################



    def simulateReads(self):

        if self.simReads == None:
            self.simReads = open(self.prefix+"_simreads.fq",'w')
            self.mutKey = open(self.prefix+"_mutation.key",'w')
            self.cntKey = open(self.prefix+"_cnt.key",'w')
                    

        for gene in self.genes:
            if not gene.seq:
                gene.getSeqFromChr(self.seq)
            gene.mutateSeq()
            self.makeReads(gene)


    def makeReads(self,gene):
  
        k=0
        exonic   = randrange(0,1000)
        intronic = randrange(0,1000)
        
        self.cntKey.write("%s EXONIC %s INTRONIC %s\n" % (gene.name,exonic,intronic))




        tCnt = len(gene.transcripts)
        readsPerTran  = int(exonic / tCnt )

        ###   MUTATE READS  - AND CALL MUTATIONS ###
        base_list = ["A","C","G","T"]
        qualStr= "".join(["B" for x in range(self.readlen)])

        for t in gene.transcripts:
            break

            tranBases = makeFeatureSeq(gene.seq,t[1][1])
            offset = len(tranBases) - self.readlen

            for i in range(readsPerTran):
                
                tIdx = randrange(offset)
                tSeq = tranBases[tIdx:tIdx+self.readlen]
                
                
                for j in range(len(tSeq)):
                    if tSeq[j] == 'x':
                        spots =  locateFeatureSpot(tIdx + j,t[1])
                        mutationInfo  =  gene.mutationKey[spots[1]]
                        if random() < mutationInfo[5]:
                            tSeq[j] = mutationInfo[4]
                            gene.mutationKey[spots[1]][7]+=1
                        else:
                            tSeq[j] = mutationInfo[3]
                            gene.mutationKey[spots[1]][6]+=1
                    if random() < 0.01:
                        tSeq[j] = choice(base_list)
            
                tmpID = "@SAMPLE_"+str(self.counter)+"_"+gene.name
                self.counter+=1

                self.simReads.write("%s\n%s\n+\n%s\n" % (tmpID,"".join(tSeq),qualStr))
        
        if len(gene.introns) > 0:
            readsPerIntron = int(intronic / len(gene.introns) )

            for intron in gene.introns:
                
                
                intTuple = intron[1][1][1]
                hgTuple  = intron[1][2][1]
                intBases = makeFeatureSeq(gene.seq,[intTuple])
               
                if len(intBases) <= self.readlen:
                    continue
                
                offset = len(intBases) - self.readlen


                for i in range(readsPerIntron):
                    intIdx = randrange(offset)
                    intSeq = intBases[intIdx: intIdx+self.readlen]

                    for j in range(len(intSeq)):

                        if intSeq[j] == 'x':
                            spots = locateFeatureSpot(intIdx + j, [[[0,len(intBases)]],[intTuple],[hgTuple]])
                            mutationInfo = gene.mutationKey[spots[1]]
                            if random() < mutationInfo[5]:
                                intSeq[j] = mutationInfo[4]
                                gene.mutationKey[spots[1]][7]+=1
                            else:
                                intSeq[j] = mutationInfo[3]
                                gene.mutationKey[spots[1]][6]+=1
                        if random() < 0.01:
                            intSeq[j] = choice(base_list)
                    
                    tmpID = "@SAMPLE_"+str(self.counter)+"_"+gene.name
                    self.counter+=1

                    self.simReads.write("%s\n%s\n+\n%s\n" % (tmpID,"".join(intSeq),qualStr))
        
        for key in sorted(gene.mutationKey):
            info = gene.mutationKey[key]
            if info != []:
                self.mutKey.write("%s %s %s %s %s %s %s %s %s %s\n" % (gene.name, info[0],gene.chr,info[1],info[2],info[3],info[4],info[5],info[6],info[7]))


##################################################################################################################################################################
       
       
