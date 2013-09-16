#!/usr/bin/env python


import sys
import os 
import difflib

from MutationFile import *
from Sequence import *
from Utilities import *
from GtfLine import *
from GtfGene import *
from random import randrange
from random import random 
from random import choice
from collections import defaultdict as dd

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix,readlen,mutationList,TYPE="ALL",printKEY=True):
       

        self.fName = open(fileHandle);  self.prefix = prefix;  self.readlen = readlen; self.PRINTKEY = printKEY
        
        self.index =0; self.genes = []; self.geneKey={}; self.seq = []
        
        self.open = True; self.FIRSTPRINT=True; self.notes = None; self.mutFile = None; self.simReads = None; self.counter = 0

        self.MUTATIONS = False
        
        tmpLine=self.fName.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr

        if mutationList != None:
            self.mutationList = True
            self.Mutations = MutationFile(mutationList)

############################################################################################################################################






    def loadGenesOnChromosome(self):
        while self.chr == self.line.chr:
            self.gene = GtfGene(self.line,self.readlen)
            while self.line.valid and self.line.geneID ==  self.gene.name:
                self.gene.addGtfLine(self.line)
                self.line = GtfLine(self.fName.readline())
            if self.gene.validOffsets():
                self.genes.append(self.gene)
             
    def startNextChromosome(self):

        if self.line.chr == 'NA':
            self.chr =  'NA'
            self.open = False
        else:
            self.chr = self.line.chr; self.genes=[]; self.seq=[]; self.geneKey={}
        

    def addFasta(self,filePath):
        c=open(filePath)
        fChr = c.readline().strip().split()[0].split(">")[1]
        #if self.chr != c.readline().strip().split(">")[1]:
        if self.chr != fChr:
            print self.chr,fChr
            
            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])

        if self.mutationList:
            while self.chr == self.Mutations.chr:
                x = self.Mutations.pos
                if self.seq[x-1] != self.Mutations.ref:
                    print "INCORRECT REF BASE"
                    print self.seq[x-1]
                    print self.Mutations.line
                    errorQuit("BAD")
                elif self.Mutations.prob > self.Mutations.minProb and self.Mutations.urn > self.Mutations.minUrn:
                    self.seq[x-1] = self.Mutations.snp
                self.Mutations.nextLine()


###########################################################################################################################################################################################################

























########################################################print###############################################################################################
################################## METHODS TO PRINT OUT ANNOTATION AND SEQUENCES ######################################################################
#######################################################################################################################################################


    def getSeqAndKey(self,gSeq,geneStart,pieces):
        mySeq = ''
        fKey = []; gKey = []; fPos = 0

        for start,end in pieces:
            fKey.append((fPos,fPos+(end-start))); gKey.append((start-geneStart,end-geneStart)); fPos += (end-start)+1
            mySeq += gSeq[start-geneStart:(end-geneStart) + 1]

        return mySeq,listToString([fKey,gKey,pieces],['-',',',' '])
                

    def printFeatures(self,gene):


        geneSeq = "".join(gene.seq)
        geneStr=gene.name+"|"+gene.hugo+"|"+gene.strand+"|"+gene.chr+"|"+gene.type+"|"+gene.status+"|"
       

        if self.printTYPE != 'EXONIC':

            for i in range(len(gene.introns)):
                seq,key = self.getSeqAndKey(geneSeq,gene.start,gene.introns[i])
                self.intronPrint.write(">%s\n%s\n" % ( geneStr+"ITRN"+str(i+1), seq))
                if self.PRINTKEY:   self.keyFile.write("%s %s\n" % (geneStr+"ITRN"+str(i+1), key))

        if self.printTYPE != 'INTRONIC':

            seq,key = self.getSeqAndKey(geneSeq,gene.start,gene.exons)
            self.exonPrint.write(">%s\n%s\n" % ( geneStr+"EXNS", seq))
            if self.PRINTKEY:   self.keyFile.write("%s %s\n" % (geneStr+"EXNS",key))

            for i in range(len(gene.knownJunks)):
                seq,key = self.getSeqAndKey(geneSeq,gene.start,gene.knownJunks[i])
                self.tranJunks.write(">%s\n%s\n" % ( geneStr+"KJXN"+str(i+1), seq))
                if self.PRINTKEY:   self.keyFile.write("%s %s\n" % (geneStr+"KJXN"+str(i+1), key))
            
            for i in range(len(gene.novelJxns)):
                seq,key = self.getSeqAndKey(geneSeq,gene.start,gene.novelJxns[i])
                self.novelJunks.write(">%s\n%s\n" % ( geneStr+"NJXN"+str(i+1), seq))
                if self.PRINTKEY:   self.keyFile.write("%s %s\n" % (geneStr+"NJXN"+str(i+1), key))
            

            seq,key = self.getSeqAndKey(geneSeq,gene.start,[(gene.start,gene.end)])
            self.genePrint.write(">%s\n%s\n" % ( geneStr+"GENE", seq))
            if self.PRINTKEY:
                key = key.split()
                key = " ".join([key[1],key[1],key[2]])
                self.keyFile.write("%s %s\n" % (geneStr+"GENE",key))

        

    def printGenesOnChromosome(self,TYPE='ALL'):
        self.printTYPE=TYPE
        if self.FIRSTPRINT:
            self.FIRSTPRINT = False
            if TYPE=='ALL' or TYPE=='EXONIC':
                self.genePrint   =   open(self.prefix+'_geneSeqs.fa','w')
                self.exonPrint   =   open(self.prefix+'_exonSeqs.fa','w')
                self.tranJunks   =   open(self.prefix+'_knownJxns.fa','w')
                self.novelJunks  =   open(self.prefix+'_novelJxns.fa','w')
            if TYPE=='ALL' or TYPE=='INTRONIC':
                self.intronPrint = open(self.prefix+'_intronSeqs.fa','w')
            if self.PRINTKEY:
                self.keyFile     = open(self.prefix+'.key','w')
        for gene in self.genes:
            if not gene.seq:
                gene.getSeqFromChr(self.seq)
            self.printFeatures(gene)















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
################################## METHODS TO SIMULATE READS FROM GTF FILE   ##########################################################################
#######################################################################################################################################################



    def simulateReads(self,mutationRate=0.01,exRate=0.01,intRate=0.01):

        if self.simReads == None:
            self.simReads = open(self.prefix+"_simreads.fq",'w')
            self.mutKey = open(self.prefix+"_mutation.key",'w')
            self.cntKey = open(self.prefix+"_cnt.key",'w')
                    

        for gene in self.genes:
            if not gene.seq:
                gene.getSeqFromChr(self.seq)
            gene.mutateSeq(exRate,intRate)
            self.makeReads(gene,mutationRate)


    def makeReads(self,gene,mutationRate):
  
        k=0
        exonic   = randrange(0,5000)
        intronic = randrange(0,1000)
        
        self.cntKey.write("%s EXONIC %s INTRONIC %s\n" % (gene.name,exonic,intronic))

        tCnt = len(gene.transcripts)
        readsPerTran  = int(exonic / tCnt )

        ###   MUTATE READS  - AND CALL MUTATIONS ###
        base_list = ["A","C","G","T"]
        qualStr= "".join(["B" for x in range(self.readlen)])


        for t in gene.transcripts:

            tranBases = makeFeatureSeq(gene.seq,t[1][1])

            if len(tranBases) <=self.readlen: continue
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
                    if random() < mutationRate:
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
                self.mutKey.write("%s %s %s %s %s %s %s %s %s %s %s\n" % (gene.name, info[1],gene.chr,info[2],gene.strand,info[0],info[3],info[4],info[5],info[6],info[7]))


##################################################################################################################################################################
       
       
