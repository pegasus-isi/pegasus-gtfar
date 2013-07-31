#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd
sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class FastqFile:
    def __init__(self,fileHandle,rlen,lowQ,avgQ,minTrim):
        self.fname = open(fileHandle)
        self.readLen ,self.lowQual, self.minAvg, self.minTrim = rlen, lowQ, avgQ, minTrim
        self.adaptors=[]
        self.multiple=10
        self.seedLen=10
        self.maxLowQuals=5
        self.adaptMatchRate=0.85
        self.outDict={}
        self.stats=False
        self.rejectFile = None
        self.rejectLabel = None
        ########## PAD THEM WITH X OR Y BUT NOT SAME DUH #############


    def recordStats(self):
        if not self.stats:
            self.stats=True
            self.allFeatures=dd(int); self.trimFeatures=dd(int); self.failFeatures=dd(int)
            self.passCnt, self.trimCnt , self.failCnt = 0,0,0
        if self.trimLoc == self.readLen:
            self.passCnt+=1
        elif self.trimLoc >= self.minTrim:
            self.trimCnt+=1
            for r in self.redFlags:
                self.allFeatures[r[1]]+=1
                self.trimFeatures[r[1]]+=1
        else:
            self.failCnt+=1
            for r in self.redFlags:
                self.allFeatures[r[1]]+=1
                self.failFeatures[r[1]]+=1

    def printStats(self):
        print "BASIC SUMMARY: PASSED/TRIMMED/REJECTED",self.passCnt,self.trimCnt,self.failCnt

        print "Filter Occurences TOTAL/RECOVERED/LOST"

        for k in sorted(self.allFeatures.keys()):
            print k,":",self.allFeatures[k],self.trimFeatures[k],self.failFeatures[k]

        
    def loadAdaptor(self,name,kind,seq):
        self.adaptors.append([name,seq,kind])
        

    def nextRead(self):
        self.readID=self.fname.readline().strip()
        self.readSeq=self.fname.readline().strip()
        self.readSign=self.fname.readline().strip()
        self.readQual=self.fname.readline().strip()
        self.phreds = [ord(q)-35 for q in self.readQual]
        

        if not self.readQual:
            self.finish=True
            return False
        else:
            self.finish=False
            
        if len(self.readSeq)!=self.readLen or len(self.readQual)!=self.readLen:
            print "ERROR"
            sys.exit()
        self.trimLoc = self.readLen;
    
        return True


    def adaptorSearch(self):
        for i in range(len(self.adaptors)):
            adaptor=self.adaptors[i]
            if adaptor[2]=="START":
                bestHit=difflib.SequenceMatcher(None,adaptor[1],self.readSeq).find_longest_match(0,self.seedLen,self.seedLen,self.readLen)
                readPiece=self.readSeq[bestHit[1]-bestHit[0]:(bestHit[1]-bestHit[0])+len(adaptor[1])]
                self.adaptors[i]=[adaptor[0],adaptor[1],adaptor[2],bestHit[1]-bestHit[0],len(readPiece),difflib.SequenceMatcher(None,readPiece,adaptor[1][0:len(readPiece)]).ratio()]
            else:
                bestHit=difflib.SequenceMatcher(None,adaptor[1],self.readSeq).find_longest_match(len(adaptor[1])-self.seedLen,len(adaptor[1]),0,self.readLen-(len(adaptor[1])-self.seedLen))
                if bestHit[0]>bestHit[1]:
                    readPiece = self.readSeq[0:len(adaptor[1])-(bestHit[0]-bestHit[1])]
                    self.adaptors[i]=[adaptor[0],adaptor[1],adaptor[2],0,len(readPiece),difflib.SequenceMatcher(None,readPiece,adaptor[1][len(adaptor[1])-len(readPiece)::]).ratio()]
                else:
                    readPiece=self.readSeq[bestHit[1]-bestHit[0]:(bestHit[1]-bestHit[0])+len(adaptor[1])]
                    self.adaptors[i]=[adaptor[0],adaptor[1],adaptor[2],bestHit[1]-bestHit[0],len(readPiece),difflib.SequenceMatcher(None,readPiece,adaptor[1][0:len(readPiece)]).ratio()]

                
    def rateQuality(self):
        self.lowQuals=[]; self.firstN=None
        for i in range(len(self.phreds)):
            if self.phreds[i]<self.lowQual and len(self.lowQuals)<self.maxLowQuals:
                self.lowQuals.append(i)
            if self.readSeq[i]=="N" and self.firstN==None:
                self.firstN=i
        self.avgQual=sum(self.phreds)/float(len(self.phreds))


    def filterPass(self):
        self.redFlags=[]
        if self.firstN!=None:
            self.redFlags.append((self.firstN,'firstN'))

        if self.avgQual<self.minAvg:
            self.redFlags.append((len(self.readSeq),'lowAvgQual'))
        
        if len(self.lowQuals)>self.maxLowQuals:
            self.redFlags.append((self.lowQuals[-1],'manyPoorQuals'))

        for adaptor in self.adaptors:
            if (adaptor[2]=="START" and adaptor[5]>0.8 and  adaptor[4]>3) or (adaptor[2]=="END" and adaptor[4]>8 and adaptor[5]>0.75):
                self.redFlags.append((adaptor[3],adaptor[0]))
      
        if len(self.redFlags)>0:
            self.redFlags.sort()
            return False
        return True
        

        
    def tryTrimFastq(self,prefix,minTrim):
        self.trimLoc=self.redFlags[0][0]
        trimPhreds=self.phreds[0:self.trimLoc+1]
        while sum(trimPhreds)/float(len(trimPhreds))<self.minAvg:
            if self.trimLoc<minTrim: break
            self.trimLoc-=1
            trimPhreds=self.phreds[0:self.trimLoc]
        if self.trimLoc>=minTrim:
            newlen=self.trimLoc-(self.trimLoc%self.multiple)
            if newlen not in self.outDict.keys():
                self.outDict[newlen]=open(prefix+'_trim_'+str(newlen)+'.fastq','w') 
            self.outDict[newlen].write("%s\n%s\n%s\n%s\n" % (self.readID,self.readSeq[0:newlen],self.readSign,self.readQual[0:newlen]))
        else:
            if self.rejectFile == None:
                self.rejectFile = open(prefix+'_reject.fastq','w')
            self.rejectFile.write('%s %s\n' % (self.readID," ".join([x[1] for x in self.redFlags])))
            self.rejectFile.write("%s\n%s\n%s\n" % (self.readSeq,self.readSign,self.readQual))
            
                
            

    def writeFastq(self,prefix):
        if self.readLen not in self.outDict.keys():
            self.outDict[self.readLen]=open(prefix+'_full_'+str(self.readLen)+'.fastq','w') 
        self.outDict[self.readLen].write("%s\n%s\n%s\n%s\n" % (self.readID,self.readSeq,self.readSign,self.readQual))

    

