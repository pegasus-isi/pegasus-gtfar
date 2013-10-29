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
    def __init__(self, fileHandle, prefix,readlen,mutationList=None,printKEY=True):
      
        try:
            self.fName = open(fileHandle)
        except TypeError:
            errorQuit("A GTF FILE IS NOT SUPPLIED")
        
        
        self.prefix = prefix;  self.readlen = readlen; self.PRINTKEY = printKEY
        
        self.index =0; self.genes = []; self.geneKey={}; self.seq = []
        
        self.open = True; self.FIRSTPRINT=True; self.notes = None; self.mutFile = None; self.simReads = None; self.counter = 0

        
        tmpLine=self.fName.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr

        if mutationList != None:
            self.mutationList = True
            self.Mutations = MutationFile(mutationList)
        else:
            self.mutationList = False

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
        try:
            c=open(filePath)
        except IOError:
            errorQuit("A valid path to chromosome fasta files was not provided")
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

    def addSpliceRecord(self,spliceRecord):
        self.spliceCands = spliceRecord.validJxns
        
########################################################################################################################################
########################################################print###############################################################################################
################################## METHODS TO PRINT OUT ANNOTATION AND SEQUENCES ######################################################################
#######################################################################################################################################################




    def featurePrint(self,features,fStr,geneStart,geneName,outPut):
        
        for F in features:
            tmpSeq= "".join([self.geneSeq[f[0]-geneStart:(f[1]-geneStart)+1] for f in F])
            if len(tmpSeq) >= self.readlen:
                outPut.write(">%s\n%s\n" % (geneName+fStr+":"+listToString(F,["|","-"]), tmpSeq))

    def printFeatures(self,gene):

        self.geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
        geneInfo = gene.name+"|"+gene.hugo+"|"+gene.type+"|"+gene.chr+"|"+gene.strand+"|"

        geneJxns = set([]);
        for g in gene.knownJxns:
            for i in range(1,len(g)): geneJxns.add(tuple((g[i-1][1],g[i][0])))

        geneJxns = sorted([g for g in geneJxns])
       
        if len(self.geneSeq) >= self.readlen:
            self.genePrint.write(">%s\n%s\n" % ( geneInfo+"GENE:"+str(gene.start)+"-"+str(gene.end), self.geneSeq))
            self.keyFile.write("%s %s \n" % (geneInfo+"GENE",listToString(gene.exons,[",","-"])))
            self.keyFile.write("%s %s \n" % (geneInfo+"JXNS",listToString(geneJxns,["|",","])))
            for start,end in gene.exons:
                if end-start > self.readlen:
                    self.exonPrint.write(">%s\n%s\n" % (geneInfo+"EXON:"+str(start)+"-"+str(end),self.geneSeq[start-gene.start:(end-gene.start)+1]))

            self.featurePrint(gene.introns,"ITRN",gene.start,geneInfo,self.intronPrint)
            self.featurePrint(gene.knownJxns,"KJXN",gene.start,geneInfo,self.exonPrint)
            self.featurePrint(gene.novelJxns,"NJXN",gene.start,geneInfo,self.exonPrint)


        

    def printGenesOnChromosome(self,TYPE='ALL'):
        if self.FIRSTPRINT:
            self.FIRSTPRINT = False
            self.genePrint   =   open(self.prefix+'_geneSeqs.fa','w')
            self.exonPrint   =   open(self.prefix+'_exonSeqs.fa','w')
            self.intronPrint = open(self.prefix+'_intronSeqs.fa','w')
            self.keyFile     = open(self.prefix+'.key','w')
        for gene in self.genes:
            self.printFeatures(gene)



    def printSpliceData(self,TYPE='ALL'):
        if self.FIRSTPRINT:
            self.FIRSTPRINT = False
           # self.genePrint   =   open(self.prefix+'_geneSeqs.fa','w')
           # self.exonPrint   =   open(self.prefix+'_exonSeqs.fa','w')
           # self.intronPrint = open(self.prefix+'_intronSeqs.fa','w')
           # self.keyFile     = open(self.prefix+'.key','w')
        for gene in self.genes:

            spliceJxns = self.spliceCands[gene.name]
            if len(spliceJxns)==0: continue

            self.geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            gStr = ">"+gene.name+"|"+gene.hugo+"|"+gene.chr+"|"+gene.strand+"|"
            
            geneInfo = ">"+gene.name+"|"+gene.hugo
            chrInfo  =  gene.chr+"|"+gene.strand
            k=1
            gene.compilePairs()
            introns=[(0,1)]+[i[1] for i in gene.introns]+[(gene.end+1,gene.end+100)]
            
            endInfo = dd(list); startInfo=dd(list)
            seqTable = dd(str)
            endTable = dd(str); startTable=dd(str)
            for s in spliceJxns:
                while k<len(introns) and introns[k][-1] < s[0]: k+=1 
                m=k
                if len(endInfo[s[0]]) == 0:
                    if introns[k-1][-1] < s[0] and s[0] <introns[k][0]:
                        if s[0] in gene.endTable:   endInfo[s[0]].append("exCANON")
                        else:                       endInfo[s[0]].append("exNOVEL")
                    else:                           endInfo[s[0]].append("intNOVEL")
                    endInfo[s[0]].append(self.geneSeq[s[0]-(gene.start+self.readlen-2):s[0]-(gene.start-1)])
                    endInfo[s[0]].append(self.geneSeq[s[0]-(gene.start-1):s[0]-(gene.start-1)+2])
                
                while m<len(introns) and introns[m][-1] < s[1]: m+=1
                if len(startInfo[s[1]])==0:
                    if introns[m-1][-1] < s[1] and s[1] < introns[m][0]:
                        if s[1] in gene.startTable: startInfo[s[1]].append("exCANON")
                        else:                       startInfo[s[1]].append("exNOVEL")
                    else:                           startInfo[s[1]].append("intNOVEL")
                    startInfo[s[1]].append(self.geneSeq[s[1]-gene.start:s[1]-gene.start+self.readlen-1])
                    startInfo[s[1]].append(self.geneSeq[s[1]-gene.start-2:s[1]-gene.start]) 
                spotStr = "-".join([str(x) for x in s])
                intStr  = "-".join([endInfo[s[0]][2],startInfo[s[1]][2]])
                mapSpots = [[s[0]-(self.readlen-2),s[0]],[s[1],s[1]+self.readlen-2]]

                if s in gene.validPairs:    jData = endInfo[s[0]][0]+"_CANON_"+startInfo[s[1]][0]
                else:                       jData = endInfo[s[0]][0]+"_NOVEL_"+startInfo[s[1]][0]


                jStr=geneInfo+"|"+jData+"|"+chrInfo+"|"+"SPLICE="+intStr+":"+listToString(mapSpots,["|","-"])
                #if s in gene.validPairs:    jStr = gStr+spotStr+"|"+intStr+"|"+endInfo[s[0]][0]+"_CANON_"+startInfo[s[1]][0]
                #else:                       jStr = gStr+spotStr+"|"+intStr+"|"+endInfo[s[0]][0]+"_NOVEL_"+startInfo[s[1]][0]
               
                sys.stdout.write("%s\n%s\n" % (jStr,"".join([endInfo[s[0]][1],startInfo[s[1]][1]])))

                
                









########################################################################################################################################################









        

       
       


#######################################################################################################################################################
################################## METHODS TO SIMULATE READS FROM GTF FILE   ##########################################################################
#######################################################################################################################################################



    def simulateReads(self,mutationRate=0.01,exRate=0.01,intRate=0.01):

        if self.simReads == None:
            self.simReads = open(self.prefix+"_simreads.fq",'w')
            self.mutKey = open(self.prefix+"_mutation.key",'w')
            self.cntKey = open(self.prefix+"_cnt.key",'w')
                    
        print len(self.genes)
        ## YOU SHOULD MUTATE THE CHR NOT THE GENE OBVI ##
        for gene in self.genes:
            self.makeReads(gene,mutationRate)
    
    def makeReads(self,gene,expressionRate):
        geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
        print geneSeq[0:5]
        print gene.name
        print gene.hugo
        print gene.start
        print len(gene.exons)
        print len(gene.introns)
        print len(gene.knownJxns)
        print len(gene.novelJxns)
        print gene.exons
        print gene.knownJxns
        print "NOVEL"
        sys.exit()


    def makeReads2(self,gene,mutationRate):
  
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
       
       
##################################### METHODS TO EVALUATE MUTATIONS  ###################################################################################
########################################################################################################################################################


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
                            





