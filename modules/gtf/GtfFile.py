#!/usr/bin/env python


import sys
import difflib
import cPickle as pickle

sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')

from tools.gtTools import seqComplement
from GtfLine import *
from GtfGene import *

##########################################################################################################################################
##################################################  MULTI-GENE CLASS START  ##############################################################
##########################################################################################################################################





class GtfFile:
    def __init__(self, fileHandle, prefix,SILENT=False):
        self.genes = []
        self.geneKeys={}
        self.seq = []
        self.index = 0
        self.notes = None
        self.readlen=100
        self.prefix = prefix
        self.fname = open(fileHandle)
        self.open = True
        tmpLine=self.fname.readline().strip()
        while tmpLine[0]=="#":
            tmpLine=self.fname.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr
            


    def loadChromosomeGenes(self):
        
        while self.chr == self.line.chr:
            
            self.gene = GtfGene(self.line,self.readlen)
            while self.line.valid and self.line.geneID ==  self.gene.name:

                self.gene.add(self.line)

                self.line = GtfLine(self.fname.readline())
            

            if self.gene.finish():
                self.geneKeys[self.gene.name]=self.gene
                self.genes.append(self.gene)
        self.index+=1
             
    
    def startNewChromosome(self):
        
        if self.line.chr == 'NA':
            self.open = False
            return 
        else:
            self.chr = self.line.chr
            self.genes = []
            self.seq=[]
            self.geneKeys={}

    def addFasta(self,filePath):
        c=open(filePath)
        if self.chr != c.readline().strip().split(">")[1]:
            print "Wrong Chromosome File Error"
            sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])

    def concatSeq(self,gSeq,restarts):
        mySeq=''
        for restart in restarts: mySeq+=gSeq[restart[0]:restart[1]+1]
        return mySeq

    def pickleOut(self):
        pickle.dump(self.exonicOutput[2],open(self.prefix+'_exonicKey.pickle','wb'))
        pickle.dump(self.intronicOutput[2],open(self.prefix+'_intronicKey.pickle','wb'))
        pickle.dump(self.catOutput[2],open(self.prefix+'_catKey.pickle','wb'))
        pickle.dump(self.fullOutput[2],open(self.prefix+'_fullKey.pickle','wb'))

    def writeOut(self,outList,name,seq,info):
        outList[0].write(">%s\n" % name)
        outList[0].write("%s\n" % seq)
        locStrs=[[] for x in range(len(info))]
        if info[0]=="GENE":
            outList[1].write("%s %s %s \n" % (name,info[1],info[2]))
            outList[2][name]=[info[1],info[2]]
        else:
            outList[2][name]=info
            for i in range(len(info)):
                for j in range(len(info[i])):
                    locStrs[i].append(str(info[i][j][0])+"-"+str(info[i][j][1]))
                locStrs[i]=','.join(locStrs[i])

            outList[1].write("%s %s\n" % (name, " ".join(locStrs)))
        
 


    def iterPrint(self,gene):

        if gene.name != "NULL":
            
            geneSeq = self.seq[gene.start-1:gene.end]
            flankSeqs = [self.seq[gene.start-gene.flankLen-1:gene.start+gene.extend] , self.seq[gene.end-gene.extend-1:gene.end+gene.flankLen]]
           
            if gene.strand=="-":
                geneSeq=seqComplement(geneSeq)
                flankSeqs= [ seqComplement(s) for s in flankSeqs[-1::-1]]

            geneSeq="".join(geneSeq).upper()
            flankSeqs = ["".join(seq).upper() for seq in flankSeqs]
            
            exCatSeq=self.concatSeq(geneSeq,gene.exons[1])
            
            ##############################################################################################################

            geneStr=gene.name+"|"+gene.hugo+"|"+gene.strand+"|"+gene.chr+"|"+gene.type+"|"+gene.status+"|"+str(len(geneSeq))+"|"+str(len(exCatSeq))+"|"
            ## PRINT FLANKS AND EXONS ##
            
            self.writeOut(self.intronicOutput,geneStr+'fpFlank',flankSeqs[0],gene.flanks[0][1])
            self.writeOut(self.intronicOutput,geneStr+'tpFlank',flankSeqs[1],gene.flanks[1][1])
            self.writeOut(self.catOutput,geneStr+'catEXON',exCatSeq,gene.exons)
                        
            ## PRINT INTRONS ##
            for i in gene.introns:
                intCatSeq=self.concatSeq(geneSeq,i[1][1])
                self.writeOut(self.intronicOutput,geneStr+"INTRON:"+str(i[0]),intCatSeq,i[1])
            
            ## PRINT TRANSCRIPTS ##
            for t in gene.transcripts:
                tranCatSeq=self.concatSeq(geneSeq,t[1][1])
                self.writeOut(self.exonicOutput,geneStr+"TRANSCRIPT:"+t[0],tranCatSeq,t[1])
                     
            ## PRINT FULL SEQS ## 
                
            self.writeOut(self.fullOutput,geneStr+"FULLSEQ",geneSeq,["GENE",str(len(geneSeq)),str(len(exCatSeq))])


        
    def cntOverLap(self,strand,ex1,ex2):
        if strand=="-":
            if ex1[1]>ex2[0] or ex2[1]>ex1[0]: return 0
            if ex1[0]>=ex2[0]:
                if ex2[1]<ex1[1]:
                    return ex2[0]-ex1[1]+1
                else:
                    return ex2[0]-ex2[1]+1
            if ex2[0]>=ex1[0]:
                if ex1[1]<ex2[1]:
                    return ex1[0]-ex2[1]+1
                else:
                    return ex1[0]-ex1[1]+1
        else:
            if ex1[1]<ex2[0] or ex2[1]<ex1[0]: return 0
            if ex1[0]<=ex2[0]:
                if ex1[1]<ex2[1]:
                    return ex1[1]-ex2[0]+1
                else:
                    return ex2[1]-ex2[0]+1
            if ex2[0]<=ex1[0]:
                if ex2[1]<ex1[1]:
                    return ex2[1]-ex1[0]+1
                else:
                    return ex1[1]-ex1[0]+1
        print "FUCK"
        sys.exit()


    def printNotes(self,pString):
        if self.notes==None:
            self.notes=open(self.prefix+'.notes','w')
        self.notes.write("%s\n" % (pString))



    def comparePrevGene(self,old,new):
        if old==None:
            return True,new
        overLap=0
        if new.strand=='-':
            old_bases = sum([X[0]-X[1]+1 for X in old.exons[2]])
            new_bases = sum([X[0]-X[1]+1 for X in new.exons[2]])
        else:
            old_bases = sum([X[1]-X[0]+1 for X in old.exons[2]])
            new_bases = sum([X[1]-X[0]+1 for X in new.exons[2]])

        for old_exon in old.exons[2]:
            for new_exon in new.exons[2]:
                overLap+=self.cntOverLap(new.strand,old_exon,new_exon)
        if overLap==0:
            return False,None
        else:
            if overLap==old_bases or overLap==new_bases:
                if old.hugo == new.hugo:
                    if old_bases>overLap:
                        self.printNotes('Warning: '+new.name+' is a subSequence of '+old.name+' and share Hugo names and will not be annotated')
                        return True,old
                    else:
                        self.printNotes('Warning: '+old.name+' is a subSequence of '+new.name+' and share Hugo name will not be annotated')
                        return True,new
                else:
                    if old_bases>overLap:
                        self.printNotes('Warning: '+new.name+' is a subSequence of '+old.name+' but remains in the annotation')
                        return False,None
                    else:
                        self.printNotes('Warning: '+old.name+' is a subSequence of '+new.name+' but remains in the annotation')
                        return False,None
            else:
                self.printNotes('Warning: '+old.name+' and '+new.name+' share '+str(overLap)+' bases ( '+str(old_bases)+' '+str(new_bases)+' )')
                return False,None

                    
    def printSeqs(self):


        if self.index == 1:

                self.exonicOutput   = [ open(self.prefix+'_exonicSeqs.fa','w'),    open(self.prefix+'_exonicKey.txt','w'),   {} ]
                self.intronicOutput = [ open(self.prefix+'_intronicSeqs.fa','w'),  open(self.prefix+'_intronicKey.txt','w'), {} ]
                self.catOutput      = [ open(self.prefix+'_concatSeqs.fa','w'),    open(self.prefix+'_concatKey.txt','w'),   {} ]
                self.fullOutput     = [ open(self.prefix+'_fullSeqs.fa','w'),      open(self.prefix+'_fullKey.txt','w'),     {} ]






        prevPos = None
        prevNeg = None

        for gene in self.genes:

            if gene.strand=="+":
                geneCompare=self.comparePrevGene(prevPos,gene)
                if geneCompare[0]==False:
                    self.iterPrint(prevPos)
                    prevPos=gene
                else:
                    prevPos=geneCompare[1]
            else:
                geneCompare=self.comparePrevGene(prevNeg,gene)
                if geneCompare[0]==False:
                    self.iterPrint(prevNeg)
                    prevNeg=gene
                else:
                    prevNeg=geneCompare[1]
        
        self.iterPrint(prevPos)
        self.iterPrint(prevNeg)





