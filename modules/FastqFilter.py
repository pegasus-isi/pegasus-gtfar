#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd


#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class FastqReads:
    def __init__(self,fileHandle,rlen,prefix): #,strictLevel):
        
        self.fileHandle,self.readLen,self.prefix = fileHandle,rlen,prefix


        self.writeTable = {self.readLen: open(self.prefix+'_full.fastq',"w"), 0: open(self.prefix+'_reject.fastq',"w") }
        self.stats = open(self.prefix+".stats","w")
        self.cnt = 0
        ### INFER FILE TYPE FROM FIRST LINE ### 
        #self.null_read = ("@DG3_NULL_READ","".join(["N" for i in range(rlen)]),"+","".join(["I" for i in range(rlen)]))
        firstLine=self.fileHandle.readline().split()
        
        if len(firstLine)==4 and firstLine[2][0]=="+":
            self.linesPerRead=1
            self.fileOpen = True
            self.nextRead = self.singleLineRead
            self.name,self.seq,self.plus,self.qual=firstLine
            if len(self.seq) < self.readLen: 
                self.dataFailure("Supplied read is shorter ("+str(len(self.seq))+") than length ("+str(self.readLen)+") parameter")
                
        else:
            self.dataFailure(self,"no multi line yet")

        

    def addTrimLengths(self,trimString):
        
        
        try:  trims = sorted([int(s) for s in trimString.split(",") if int(s) < self.readLen and int(s)>=50])
        except ValueError:  self.dataFailure("Trim lengths supplied contain non-integer value")
        if len(trims)>0:
            if trims[-1] >= self.readLen: self.dataFailure("Trim lengths supplied contain value equal to or exceeding read length") 
            for t in trims:    self.writeTable[t] = open(self.prefix+"_"+str(t)+".fastq","w")
        
        
        self.trimHash  = {}
        trims.append(self.readLen)
        for i in range(len(trims)-1):   self.trimHash.update([(j,trims[i]) for j in range(trims[i],trims[i+1])])
       





    def addAdaptors(self,tpAdaptor,fpAdaptor,strictness):
        self.tpAdaptors,self.fpAdaptors = [s for s in tpAdaptor.split(',')],[s for s in fpAdaptor.split(',')]   
        self.tpLens = [len(x) for x in self.tpAdaptors]; self.fpLens = [len(x) for x in self.fpAdaptors]
        if len(self.tpAdaptors) != 2 or len(self.fpAdaptors) != 2:                          self.dataFailure("Adaptors must include two pieces w/ min lengths of 12bp")
        elif min([len(a) for a in self.tpAdaptors]+[len(b) for b in self.fpAdaptors]) < 12: self.dataFailure("Adaptors must include two pieces w/ min lengths of 12bp")
        elif strictness <1 or strictness > 5:                                               self.dataFailure("Supplied strictness must be between 1 and 5")
        else:
            self.minMatch = 7 + strictness 
            self.minScr   = [0.8, 0.825, 0.80, 0.875, 0.90, 0.925, 0.95][strictness]
        
            self.fpEndLocs   = [0 for i in range(self.readLen+1)]
            self.tpStartLocs = [0 for i in range(self.readLen+1)]
            self.fpCnt, self.tpCnt, self.doubleCnt, self.fpTrim, self.tpTrim, self.doubleTrim = 0,0,0,0,0,0
    

    def singleLineRead(self):
        try: self.name,self.seq,self.plus,self.qual=self.fileHandle.readline().split()
        except ValueError:
            self.fileOpen = False 
        self.cnt +=1
        


    def queryScr(self,query,qLen):
         
        qFind=difflib.SequenceMatcher(None,query,self.seq).find_longest_match(0,qLen,0,self.readLen)
        if qFind[2] >= self.minMatch: 
            seqIdxs   = max(0,qFind[1]-qFind[0]),min(self.readLen,qFind[1]+(qLen-qFind[0]))
            qIdxs = max(0,qFind[0]-qFind[1]),min(qLen,self.readLen-qFind[1]+1)
            if difflib.SequenceMatcher(None,self.seq[seqIdxs[0]:seqIdxs[1]],query[qIdxs[0]:qIdxs[1]]).ratio() >= self.minScr: 
                self.adaptorLocs = seqIdxs
                return True
        return False        
        




    def adaptorFilter(self):
    
        ## STARTING WTIH THREE PRIME HEAD TO TAIL ##
        ### FIRST SEARCH FOR FIVE PRIME ADAPTORS --- START WITH THE END ---- ####
        
        if self.queryScr(self.fpAdaptors[1],self.fpLens[1]):    fp = self.adaptorLocs[1]
        elif self.queryScr(self.fpAdaptors[0],self.fpLens[0]):  fp = self.readLen
        else:                                                   fp = 0

        if self.queryScr(self.tpAdaptors[0],self.tpLens[0]):    tp = self.adaptorLocs[0]
        elif self.queryScr(self.tpAdaptors[1],self.tpLens[1]):  tp = 0 
        else:                                                   tp = self.readLen    

            
        if fp == 0 and tp == self.readLen: self.writeTable[self.readLen].write("%s\n%s\n%s\n%s\n" % (self.name, self.seq[0:self.readLen], self.plus, self.qual[0:self.readLen]))  
        
        elif fp == 0 and tp < self.readLen:
            self.tpCnt+=1
            self.tpStartLocs[tp]+=1
            try:
                self.writeTable[self.trimHash[tp]].write("%s %s\n%s\n%s\n%s\n" % (self.name,"TPA,"+str(tp),self.seq[0:tp],self.plus,self.qual[0:tp]))
                self.tpTrim+=1             
            except KeyError: 
                self.writeTable[0].write("%s %s %s %s\n%s\n%s\n%s\n" % (self.name,"TPA",tp,self.seq[tp:self.readLen],self.seq,self.plus,self.qual))
        elif fp > 0 and tp == self.readLen:
            self.fpCnt += 1
            self.fpEndLocs[fp]+=1 
            try:
                self.writeTable[self.trimHash[self.readLen-fp]].write("%s %s\n%s\n%s\n%s\n" % (self.name,"FPA,"+str(fp),self.seq[fp:self.readLen],self.plus,self.qual[fp:self.readLen]))
                self.fpTrim+=1  
            except KeyError: self.writeTable[0].write("%s %s %s %s\n%s\n%s\n%s\n" % (self.name,"FPA",fp,self.seq[0:fp],self.seq,self.plus,self.qual))
        else:
            self.doubleCnt+=1
            self.tpStartLocs[tp]+=1; self.fpEndLocs[fp]+=1
            splits = [self.seq[0:fp],self.seq[fp:tp],self.seq[tp:self.readLen]]
            try:  
                self.writeTable[self.trimHash[len(splits[1])]].write("%s %s\n%s\n%s\n%s\n" % (self.name,"AA,"+str(fp)+","+str(tp),splits[1],self.plus,self.qual[fp:tp]))
                self.doubleTrim+=1
            except KeyError: self.writeTable[0].write("%s %s %s %s %s %s\n%s\n%s\n%s\n" % (self.name,"AA",fp,tp,splits[0],splits[2],self.seq,self.plus,self.qual[fp:self.readLen]))

            
            

            
    def dataFailure(self,msg="FOO"):
        sys.stderr.write(msg+"\n")
        sys.exit(2)



    def printSummaryStats(self):
        self.stats.write( "-info1-TOTAL_READS %s\n" % self.cnt)
        self.stats.write( "-info2-PASSED %d\n"  %   (self.cnt - (self.tpCnt+self.fpCnt+self.doubleCnt)))
        self.stats.write( "-info3-TRIMMED %s\n" %    (self.tpTrim + self.fpTrim + self.doubleTrim))
        self.stats.write( "-info4:THREE_PRIME_ADAPTOR %s\n" % self.tpCnt)
        self.stats.write( "-info5:THREE_PRIME_TRIM %s\n" % self.tpTrim) 
        self.stats.write( "-info6:FIVE_PRIME_ADAPTOR %s\n" % self.fpCnt)
        self.stats.write( "-info7:FIVE_PRIME_TRIM %s\n"  % self.fpTrim)
        self.stats.write( "-info8:DOUBLE_ADAPTOR %s\n" % self.doubleCnt)
        self.stats.write( "-info9:DOUBLE_TRIM %s\n" % self.doubleTrim)
        for i in range(0,self.readLen):
            if i<9: self.stats.write( "-locs:fpEnd/tpStart:00"+str(i+1)+" %s %s\n" % (self.fpEndLocs[i],self.tpStartLocs[i]))
            elif i<99: self.stats.write( "-locs:fpEnd/tpStart:0"+str(i+1)+" %s %s\n" % (self.fpEndLocs[i],self.tpStartLocs[i]))
            else: self.stats.write( "-locs:fpEnd/tpStart:"+str(i+1)+" %s %s\n" % (self.fpEndLocs[i],self.tpStartLocs[i]))
            
