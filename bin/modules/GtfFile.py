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
from random import randrange
from random import random 
from random import choice
from collections import defaultdict as dd

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix=None,readlen=None,filterType=None,findCands=True):
        try:    self.fName = open(fileHandle)
        except TypeError:   errorQuit("A GTF FILE IS NOT SUPPLIED")
            
        self.prefix, self.readLen,self.findCands = prefix, readlen,findCands
     
        self.genes, self.seq = [], [] 
        self.open = True
        self.minLen,self.maxLen = 35,400 
        self.SNPCANDS = False
        
        if self.prefix != None:
            self.featureFile =   open(self.prefix+'_features.fa','w')
            self.chrFile =   open(self.prefix+'_chrs.fa','w')
            self.headerFile  =   open(self.prefix+"_headers.txt","w")
            self.headerFile.write("%s\n" % "\t".join(["@HD","VN:0.1.5c","SO:queryname"]))
            self.headerFile.write("%s\n" % "\t".join(["@SQ","SN:chrR","LN:10000"]))
            if self.findCands:
                self.geneFile    =   open(self.prefix+'_geneSeqs.fa','w')
                self.candFile    =   open(self.prefix+'_jxnCands.fa','w')
       
         
        
        
        if filterType == "HUMAN":
            for name,seq in GtfFilters(filterType).seqs: self.featureFile.write("%s:0-%s\n%s\n" % (name,len(seq),seq))

        tmpLine=self.fName.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr



############################################################################################################################################
######################################## Annotation Methods ################################################################################
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
            if len(filePath.split("/")) > 1:
                failChr = filePath.split("/")[-1]
            else: failChr = filePath
            errorQuit("Error: "+failChr+" is referenced by the suppled gtf-file but not found in the supplied chromosome path: "+filePath)
        fChr = c.readline().strip().split()[0].split(">")[1]
        #if self.chr != c.readline().strip().split(">")[1]:
        if self.chr != fChr:
            print self.chr,fChr
            
            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])

    
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
        
        
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

                
############################################################################################################################################
####################################################### CANDIDATE CODE  ####################################################################
############################################################################################################################################
                
    def addCandidates(self,candidate_file,cand_type):
        if cand_type=="SNPS":
            self.SNPCANDS = True
            self.snpcands = dd(lambda: dd(list))
            for line in open(candidate_file):
                line=line.split()
                self.snpcands[line[0]][int(line[1])]=[line[2::]]
                


    def printCandidatesOnChromosome(self,TYPE='ALL'):
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
        
        
############################################################################################################################################





############################################################################################################################################
####################################################### SIMULATION CODE ####################################################################
############################################################################################################################################

    def add_simulation_parameters(self,simulation_parameter_file):
        for line in open(simulation_parameter_file):
            line = line.split()
            if line[0] == 'cdna_length':         self.cdna_length  = int(line[-1])
            if line[0] == 'sonication_mean':     self.sonication_mean = int(line[-1])
            if line[0] == 'tail_adaptor'    :    self.tail_adaptor = line[-1]
            if line[0] == 'head_adaptor'    :    self.head_adaptor = line[-1]
            if line[0] == 'output_file_prefix':  self.prefix = line[-1]
            if line[0] == 'read_length'    :     self.readLen = int(line[-1])
            if line[0] == 'protocol'    :     self.protocol  = line[-1]
            if line[0] == 'size_selection_min'    :     self.size_selection_min  = int(line[-1])
            if line[0] == 'linear_amplification_values': self.linear_amplification_values = [int(x) for x in line[2::]]
            if line[0] == 'exponential_amplification_values': self.exponential_amplification_values = [int(x) for xi in line[2::]]
    



    def add_gene_key(self,gene_key_file):
        self.gene_key=dd(lambda: (1,0))
        for line in open(gene_key_file):
            line = line.split()
            self.gene_key[line[0]]=(int(line[-2]),int(line[-1]))











    def simulateGeneReadsOnChromosome(self,protocol="HUGO"):
        
        self.read_noise   = "".join(['A' for i in range(self.readLen)])
        self.qual = "".join(["I" for i in range(self.readLen)])
        self.readOutPut = sys.stdout
        
        for gene in self.genes:
            read_num = 0
            vals = self.gene_key[gene.hugo]
            gene.fragmentRNA(self.cdna_length,vals[0],vals[1])
            if self.protocol == "HUGO" or self.protocol == 'hugo':
                gene.linearAmplifyFrags(self.linear_amplification_values)
                gene.sonicateFrags(self.sonication_mean)
          
            if self.protocol == "OLEG" or self.protocol == 'oleg':
                gene.sonicateFrags(self.sonication_mean)
                gene.exponentialAmplifyFragments(self.exponential_amplification_values)
                 
           
            
            geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            for frag in gene.rnaFragments:
                if frag[-1] < self.size_selection_min: continue 
                read_num+=1
                readID = ["@"+str(read_num),gene.chr,gene.name,frag[1],self.protocol]
                readSeq = "".join([geneSeq[x[0]-gene.start : (x[1]-gene.start)+1] for x in frag[2]])[0:self.readLen]
                read_coords= []
                read_dist,k   = 0,0
                while k < len(frag[2]):
                    if read_dist + (frag[2][k][1] - frag[2][k][0] ) + 1 <= len(readSeq):
                        read_coords.append(frag[2][k])
                        read_dist += (read_coords[-1][1] - read_coords[-1][0] ) + 1
                    else:
                        read_coords.append((frag[2][k][0],frag[2][k][0] + (len(readSeq) - read_dist)))
                        break
                    k+=1
                coord_strs = "G".join([str(x[0])+'-'+str(x[1]) for x in read_coords])
                readID.append(coord_strs)
                if len(readSeq) == self.readLen:    readID.append("N")
                else: 
                    readID.append("Y")
                    if len(readSeq) < self.readLen: readSeq += self.tail_adaptor[0 : (self.readLen - len(readSeq))  ] 
                    if len(readSeq) < self.readLen: readSeq = self.head_adaptor[max(0,len(self.head_adaptor)-(self.readLen-len(readSeq))) : len(self.head_adaptor) ]  + readSeq
                    if len(readSeq) < self.readLen: readSeq += self.read_noise[0: (self.readLen - len(readSeq)) ]
                
                for i in range(frag[0]):
                    #readID.append(str(i)) 
                    tmpID = "_".join(readID + [str(i)])
                    #readID_tmp = "_".join([j for j in readID])       
                    self.readOutPut.write("%s\n%s\n%s\n%s\n" % (tmpID,readSeq,'+',self.qual))

        




        

       
       


