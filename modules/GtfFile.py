#!/usr/bin/env python


import sys
import os 
import difflib

#from MutationFile import *
#from Sequence import *
from ToolSet import errorQuit
from ToolSet import listToString
from GtfLine import *
from GtfGene import *
from GtfFilters import * 
#from random import randrange
#from random import random 
#from random import choice
from collections import defaultdict as dd

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix,readlen,filterType=None,findCands=True):
      
        try:    self.fName = open(fileHandle)
        except TypeError:   errorQuit("A GTF FILE IS NOT SUPPLIED")
       
        self.prefix, self.readLen,self.findCands = prefix, readlen,findCands
        self.genes, self.seq = [], [] 
        self.open = True
        self.minLen,self.maxLen = 35,400 
        
        self.featureFile =   open(self.prefix+'_features.fa','w')
        self.chrFile =   open(self.prefix+'_chrs.fa','w')
        self.headerFile  =   open(self.prefix+"_headers.txt","w")
        self.headerFile.write("%s\n" % "\t".join(["@HD","VN:0.1.5c","SO:queryname"]))
        self.headerFile.write("%s\n" % "\t".join(["@SQ","SN:chrR","LN:10000"]))
        
        
        
        if self.findCands:
            self.geneFile    =   open(self.prefix+'_geneSeqs.fa','w')
            self.candFile    =   open(self.prefix+'_jxnCands.fa','w')
        #self.interGenic  =   open(self.prefix+'_genome.fa','w')
        #self.novelJxns   =   open(self.prefix+'_novelJxns.fa','w')
        #self.jnxCands   =   open(self.prefix+'_novelJxns.fa','w')
        if filterType == "HUMAN":
            for name,seq in GtfFilters(filterType).seqs: self.featureFile.write("%s:0-%s\n%s\n" % (name,len(seq),seq))
       # for name,seq in self.filters.seqs: self.featureFile.write("%s:0-%s\n%s\n" % (name,len(seq),seq))

        tmpLine=self.fName.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr


############################################################################################################################################

    def loadGenesOnChromosome(self):
        while self.chr == self.line.chr:
            self.gene = GtfGene(self.line,self.readLen)
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
        try:
            c=open(filePath)
        except IOError:
            print filePath
            errorQuit("A valid path to chromosome fasta files was not provided")
        fChr = c.readline().strip().split()[0].split(">")[1]
        #if self.chr != c.readline().strip().split(">")[1]:
        if self.chr != fChr:
            print self.chr,fChr
            
            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])


        
########################################################################################################################################
########################################################print###############################################################################################
################################## METHODS TO PRINT OUT ANNOTATION AND SEQUENCES ######################################################################
#######################################################################################################################################################




    
    def printGenesOnChromosome(self,TYPE='ALL'):
        geneTuples=[]
        for gene in self.genes:
            geneTuples.append((gene.start-1,gene.end))
            geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            geneInfo = gene.name+"|"+gene.hugo+"|"+gene.type+"|"+gene.chr+"|"+gene.strand+"|"
            if gene.end - gene.start > self.minLen: 
                if self.findCands:
                    self.geneFile.write(">%s\n%s\n" % ( geneInfo+"GENE:"+str(gene.start)+"-"+str(gene.end), geneSeq)) 
                    for N in gene.novelJxns: self.candFile.write(">%s\n%s\n" % (geneInfo+"NJXN:"+listToString(N,["|","-"]),"".join([geneSeq[n[0]-gene.start:(n[1]-gene.start)+1] for n in N])))
            
            
            
            for start,end in gene.exons:    self.featureFile.write(">%s\n%s\n" % (geneInfo+"EXON:"+str(start)+"-"+str(end),geneSeq[start-gene.start:(end-gene.start)+1]))
            for I in gene.introns:   self.featureFile.write(">%s\n%s\n" % (geneInfo+"INTRON:"+listToString(I,["|","-"]),"".join([geneSeq[i[0]-gene.start:(i[1]-gene.start)+1] for i in I])))
            for J in gene.knownJxns: self.featureFile.write(">%s\n%s\n" % (geneInfo+"KJXN:"+listToString(J,["|","-"]),"".join([geneSeq[j[0]-gene.start:(j[1]-gene.start)+1] for j in J])))
            
        
        for g in geneTuples:
            if (g[1]-g[0]) > (self.readLen+1)*2:
                self.seq[g[0]+self.readLen:g[1]-self.readLen]=["N" for i in range((g[1]-self.readLen)-(g[0]+self.readLen))]
        
        chrLen=len(self.seq)


        self.headerFile.write("%s\n" % "\t".join(["@SQ","SN:"+self.chr,"LN:"+str(chrLen)]))
        k=0;BUFFER=200
        self.chrFile.write("%s\n" % (">"+self.chr))
        while True:
            self.chrFile.write("%s\n" % "".join([b for b in self.seq[BUFFER*k:BUFFER*(k+1)]]))
            k+=1
            if k*BUFFER > chrLen: break
        
        
        #self.chrFile.write("%s\n%s\n" % (">"+self.chr,self.seq))

                
                









########################################################################################################################################################









        

       
       


