#!/usr/bin/env python


import sys
from MapLine import *
from MapRead import *
from Sequence import *
from Utilities import *
from Bio.Seq import reverse_complement
from collections import defaultdict as dd
from math import fabs

##########################################################################################################################################
#####################################################  MAPDATA CLASS START      ##########################################################
########print#############################################################################################################################

class MapData:

    def __init__(self,prefix,fName):
 


        self.prefix, self.refType,self.mapFile = prefix,"EXONIC",fName
        
        self.spliceExp, self.chrCnt, self.groupCnts, self.subCnts = dd(int), dd(int), dd(int), dd(int)
       
        self.senseCnt, self.mitoCnt, self.ambigCnt , self.multiAnnoCnt, self.repetitiveCnt, self.hgUniqCnt, self.geneUniqCnt, self.readCnt, self.spliceCnt,self.uniqCnt = 0,0,0,0,0,0,0,0,0,0
        
        self.geneCnts = dd(lambda: [0,0]); self.multiGeneCnts = dd(lambda: [0,0]); self.spliceCnts = dd(int); self.novelSpliceCnts = dd(int)
        
        self.geneKey, self.biasTable  = dd(lambda: ["MULTI",[0,0,0,0]]), dd(lambda:[0 for i in xrange(100)])

        self.switchStrand = {'0':'+','16':'-','+':'0','-':'16'}



        if self.refType != "INTERGENIC":
                
            
            self.samOut      = open(self.prefix+'.vis','w')
            self.geneOut     = open(self.prefix+'.genecnts','w') 
            self.spliceOut   = open(self.prefix+'.splicecnts','w')
            self.process     = self.processGeneRead
            self.printResult = self.printGeneResult
        else:
            self.process = self.processGenomeRead
            self.printResult = self.printGenomeResult

##########################################################################################################################################
#####################################################  INIT FUNCTIONS  ref         ##########################################################
#########################################################################################################################################
 


    def determineSplices(self,mapRead):
        if not mapRead.hgUniq:  return None
        



    def processGeneRead(self,mapRead):
        self.readCnt +=1
        self.subCnts[mapRead.subs]+=1
        if mapRead.sense:   self.senseCnt +=1
        if mapRead.mito:    self.mitoCnt +=1
        if mapRead.spliced:
            self.spliceCnt+=1
            self.spliceCnts[mapRead.spliceSite]+=1

        
        ###################################################################################################3
        if mapRead.hgUniq:
            self.hgUniqCnt +=1
            self.samOut.write("%s\n" % self.samString(mapRead.chr,mapRead.strand,mapRead.locs,mapRead,"UNIQUE"))


        elif mapRead.repetitive:
            self.repetitiveCnt +=1
            for m in mapRead.multiLocs:
                self.samOut.write("%s\n" % self.samString(m.chr,m.strand,m.locs,mapRead,"REPETITIVE"))

        else:
            ## ASSUMING WE ARE AMBIGUOUS HERE ##
            if not mapRead.ambiguous:
                print "WTF"
            self.ambigCnt +=1

        ##########################################################################################################

        if mapRead.geneUniq:
            self.geneUniqCnt+=1
            self.geneCnts[mapRead.geneID+"|"+mapRead.geneAltID][0]+=1
            self.groupCnts[mapRead.geneGroup]+=1
        else:
            if mapRead.hgUniq:
                self.multiAnnoCnt+=1
                self.geneCnts[",".join([g.geneID+"|"+g.geneAltID for g in mapRead.multiGenes])][0]+=1
            if mapRead.repetitive: 
                for g in mapRead.multiGenes:
                    self.geneCnts[g.geneID+"|"+g.geneAltID][1]+=1
        






    def processGenomeRead(self,mapRead):
        if mapRead.invalid:
            self.invalidCnt +=1 
            return
        self.readCnt +=1
        self.senseCnt += 1

        if mapRead.hgUniq:
            mapRead.hgJxnSearch()
            if mapRead.spliced:
                self.spliceCnt+=1
                for s in mapRead.spliceJxns: self.spliceExp[s]+=1
            for read in mapRead.splitReads:
                sys.stdout.write("%s\n" % " ".join(read))

            #print "YO",mapRead.hgLoc,mapRead.geneLoc,mapRead.read







    def samString(self,chr,strand,myPos,mapRead,TAG):
        
        cigar = "".join([str(myPos[i]-myPos[i-1]+1)+"M" if i%2==1 else str(myPos[i]-myPos[i-1]-1)+"N" for i in xrange(1,len(myPos))])
        return "\t".join([mapRead.name,self.switchStrand[strand],chr,str(myPos[0]),'255',cigar,'*','0','0',mapRead.seq,mapRead.qual,'NM:i:'+str(mapRead.subs),"GT:i:"+TAG])

###############################################################################################################################################################33







    def printGeneResult(self):

        for gene in self.geneCnts:
            self.geneOut.write("%s %s %s\n" % (gene,self.geneCnts[gene][0],self.geneCnts[gene][1]))
        for splice in self.spliceCnts:
            self.spliceOut.write("%s %s\n" % (splice,self.spliceCnts[splice]))

        mapSubs = [str(x[0])+": "+str(x[1]) for x in sorted([(k,self.subCnts[k]) for k in self.subCnts])]
        
        self.statsOut = open(self.prefix+'.stats','w')
        
        self.statsOut.write("QC-Statistics: %s \n" % self.mapFile)
        self.statsOut.write("Total-Reads[total,uniq/repetitive/ambiguous] %s %s %s %s\n" % (self.readCnt, self.hgUniqCnt,self.repetitiveCnt,self.ambigCnt))
        self.statsOut.write("Total-Reads[Mismatches] %s\n" % (" ".join(mapSubs)))
        self.statsOut.write("Uniq-Reads[spliced/unspliced] %s %s \n" % (self.spliceCnt ,self.hgUniqCnt-self.spliceCnt))
        self.statsOut.write("MitoChondrial-Reads %s \n" % self.mitoCnt)
       
        
        for t in self.groupCnts:  self.statsOut.write("GROUP: %s-Reads %s\n" % (t.upper(),self.groupCnts[t]))

         
            

    def printGenomeResult(self):


        if self.spliceCnt > 0:
            self.spliceOut = open(self.prefix+'_splice.cnts','w')
            for s in self.spliceExp:
                self.spliceOut.write("%s %s %s\n" % (s,self.refType,self.spliceExp[s]))

        self.statsOut = open(self.prefix+'.stats','w')
        self.statsOut.write("QC-Statistics: %s \n" % self.mapFile)
        self.statsOut.write("Total-Reads[valid/invalid] %s %s \n" % (self.readCnt, self.invalidCnt))
        self.statsOut.write("Valid-Reads[sense/antisense] %s %s \n" % (self.senseCnt , self.readCnt-self.senseCnt))
        self.statsOut.write("Valid-Reads[spliced/unspliced] %s %s \n" % (self.spliceCnt , self.readCnt-self.spliceCnt))
        self.statsOut.write("MitoChondrial-Reads %s \n" % self.chrCnt["chrM"])
        

         
            


##############################################################################################################
####################################   PRINTING DATA #########################################################
##############################################################################################################

    






































































                
            


    def close(self):

        self.mapHandle.close()


        if self.samOut: self.samOut.close()
        if self.mapOut: self.mapOut.close()
        if self.geneMaps > 0: self.geneOut.close()
        if self.spliceMaps > 0: self.spliceOut.close()
        #if self.novelMaps > 0: self.novelOut.close()






