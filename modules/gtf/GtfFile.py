#!/usr/bin/env python


import sys
import os 
import difflib
import cPickle as pickle


from ..tools.gtTools import *
from GtfLine import *
from GtfGene import *

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix,readlen,printKEY):
        
        self.fName = open(fileHandle);  self.prefix = prefix;  self.readlen = readlen; self.PRINTKEY = printKEY
        
        self.index =0; self.genes = []; self.geneKey={}; self.seq = []
        
        self.open = True; self.FIRSTPRINT=True; self.notes = None
        
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
        self.index+=1
             
    
    def startNextChromosome(self):
        
        if self.line.chr == 'NA':
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




    def printAnnotation(self,TYPE='ALL'):
        self.printTYPE=TYPE

        if self.FIRSTPRINT:
            self.FIRSTPRINT = False
            if TYPE=='ALL' or TYPE=='EXONIC':
                self.genePrint   = open(self.prefix+'_geneSeqs.fa','w')
                self.exonPrint   = open(self.prefix+'_exonSeqs.fa','w')
                self.catsOnly   = open(self.prefix+'_catsOnly.fa','w')
            if TYPE=='ALL' or TYPE=='INTRONIC':
                #self.genePrint   = open(self.prefix+'_geneSeqs.fa','w')
                self.intronPrint = open(self.prefix+'_intronSeqs.fa','w')
            
            
            if self.PRINTKEY:
                self.keyFile     = open(self.prefix+'.key','w')

        for gene in self.genes:
            if not gene.seq:
                gene.getSeqFromChr(self.seq,0)
            gene.forceIndex(0)
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
            if self.PRINTKEY:
                self.keyFile.write("%s %s\n" % (catString,self.junctionLocations(gene.exons)))
            for t in gene.transcripts:
                tranSeq= self.concatSeq(geneSeq,t[1][1]); tCnt+=1
                tranStr=geneStr+"TRAN"+str(tCnt)+":"+str(len(tranSeq))+"|"
                self.exonPrint.write(">%s\n%s\n" % (tranStr+t[0],tranSeq))
                if self.PRINTKEY: self.keyFile.write("%s %s\n" % (tranStr,self.junctionLocations(t[1])))
        
            self.genePrint.write(">%s\n%s\n" % (catString+"FULLSEQ",geneSeq))
        
        


       
       








