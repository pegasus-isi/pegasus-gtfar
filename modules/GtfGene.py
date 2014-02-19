#!/usr/bin/env python


import sys
from random import randrange
from random import choice
from collections import defaultdict as dd

##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class GtfGene:
    def __init__(self,line,readlen=100,flankLen=200,minOverlap=5):
        if line==None:
            self.name = None
        elif line=="NULL":
            self.start=0; self.end=0; self.name="NULL"
        elif line.type != 'gene':
            print line.type,"YO"
            print line.geneID, line.geneType,line.hugoName
            print "BAD INIT"
            sys.exit()
        else:
            self.name, self.chr, self.start, self.end, self.strand, self.type, self.status,self.hugo = line.geneID, line.chr, line.start,line.end,line.strand,line.geneType,line.geneStatus,line.hugoName 
            self.readlen = readlen; self.length = self.end-self.start+1
            
            self.transcripts = []
            self.EXONS = set([])
            self.seq = []
            self.maxJump = 10
            self.mutationKey = None
            
###############################print###########################################################################################################

#1)  ADDING ANOTHER LINE FROM A GTF FILE [ NOTICE YOU MUST START WITH GENE, THEN TRANSCRIPT OR EXON OK  


    def addGtfLine(self, line):
        if line.type=='transcript':
            self.transcripts.append([line.tranID,[]])
        elif line.type == 'exon':
            if line.start < self.start:
                self.GENOME_RELATIVE = False
            self.transcripts[-1][1].append((line.start,line.end))
            self.EXONS.add((line.start,line.end))
        else:
            return
  

###########################################################################################################################################



    def checkOffsets(self):

        self.exons = sorted(list(self.EXONS))
        
        if self.exons[0][0] < self.start:
            for i in range(len(self.exons)):
                self.exons[i] = (self.exons[i][0] + self.start -1, self.exons[i][1] + self.start -1)
            for i in range(len(self.transcripts)):
                self.transcripts[i][1].sort()
                for j in range(len(self.transcripts[i][1])):
                    self.transcripts[i][1][j] = (self.transcripts[i][1][j][0] + self.start -1, self.transcripts[i][1][j][1] + self.start -1)
        if self.exons[-1][1] > self.end:
            self.end = self.exons[-1][1]




    def findMajorExons(self):
        majorExons, majorSplice , myDists, tmpJxns = [],[], [], []
        k=0; OPEN_FRAME = True
        exStart=self.exons[0][0]; exEnd=self.exons[0][1]
        for i in range(1,len(self.exons)):
            if self.exons[i][0] <= exEnd:
                if self.exons[i][1] > exEnd: exEnd = self.exons[i][1]
            else:
                majorExons.append((exStart,exEnd))
                exStart=self.exons[i][0]; exEnd=self.exons[i][1]
        
        majorExons.append((exStart,exEnd))
        self.exons = sorted(majorExons)




#######################################################################################################################################################


    

    def findKnownJunctions(self):
        tmpJxns=set([]); maxlen = self.readlen - 1 
        tmpShort = set([]); tmpLong = set([]); tmpMini=set([])
        for t in self.transcripts:
            if len(t[1])<2: continue
            k=1;
            tJxns = sorted(t[1])
            tExons = sorted(t[1])
            tDists = [(x[1]-x[0])+1 for x in tExons]
            if tDists[0] >= maxlen:
                PREV_MAJOR=True
                MAJOR=True; MINOR=False
                mySplice = [(tExons[0][1]-(maxlen)+1,tExons[0][1])]
                myDists  = [maxlen]
            else:
                PREV_MAJOR = False
                MAJOR=False; MINOR=True
                mySplice = [(tExons[0][0],tExons[0][1])]
                myDists  = [tDists[0]]
            for i in range(1,len(tJxns)):

                if tDists[i] >= maxlen:
                    tmpShort.add((mySplice[-1],(tExons[i][0],tExons[i][0]+maxlen-1)))
                    if len(mySplice) > 1:
                        tmpLong.add(tuple( (mySplice+[(tExons[i][0],tExons[i][0]+(maxlen-myDists[-1]-1))]) ))
                    mySplice = [(tExons[i][1]-(maxlen)+1,tExons[i][1])]
                    myDists  = [maxlen]
                    PREV_MAJOR = True
                else:
                    if PREV_MAJOR:
                        tmpShort.add(tuple((mySplice + [tExons[i]])))
                        mySplice = [(mySplice[-1][1]-(maxlen-tDists[i]-1),mySplice[-1][1]),tExons[i]]
                        ## IN THIS CASE I DONT WANT TO ADD ## - BUT YOU WONT CAUSE ITS NOT ----GREATER---- THAN MAXLEN
                        myDists = [x[1]-x[0]+1 for x in mySplice]
                        PREV_MAJOR = False
                    else:
                        if sum(myDists[1::])+tDists[i] >= maxlen:
                            if len(mySplice)<2:
                                print "WTF"
                                sys.exit()
                            #if tExons[i][0] == 1256376: print "YO",mySplice ,tExons[i],i+1,len(tJxns)
                            if i+1 < len(tJxns):
                                tmpLong.add(tuple ( mySplice + [(tExons[i][0],tExons[i][0]+(maxlen-sum(myDists[1::])-1))]) )
                                newSplice = [mySplice[j] for j in range(len(mySplice)) if sum(myDists[j+1::])+tDists[i] < maxlen]+[tExons[i]]
                                newDists  = [x[1]-x[0]+1 for x in newSplice]
                                mySplice = [(newSplice[0][1]-(maxlen-sum(newDists[1::])-1),newSplice[0][1])] +newSplice[1::]
                                myDists = [x[1]-x[0]+1 for x in mySplice]
                           
                            else:
                                tmpLong.add(tuple( mySplice + [tExons[i]]))
                                break
                        else:
                            mySplice = mySplice+[tExons[i]]
                            myDists = [x[1]-x[0]+1 for x in mySplice]
                            if i+1 == len(tJxns):
                                if len(mySplice) == 2: tmpShort.add(tuple(mySplice))
                                elif len(mySplice) >2: tmpLong.add(tuple(mySplice))
                                break

        passShort = self.uniquifyEnds(tmpShort)
        passLong = self.uniquifyEnds(tmpLong)
        self.knownJxns = sorted(passShort + passLong)
















    def findNewJunctions(self):
        self.novelJxns=[]; novelSet = set([]); maxlen = self.readlen - 1
        tmpEnds,tmpStarts,jPoints = [],[],[]
        prevSites = set([])
        tmpEnds = {}; tmpStarts = {}
        jxnEnds = set([]); jxnStarts = set([]); obsSites = []
        for p in self.knownJxns:
            obsSites += [(p[i-1][1],p[i][0]) for i in range(1,len(p))]
            if len(p)==2:
                jxnEnds.add(p[0:1])
                jxnStarts.add(p[1:2]) 
            elif len(p)==3:
                jxnEnds.add(p[0:2])
                jxnStarts.add(p[0:2])
            elif len(p)>3:
                for i in range(2,len(p)-1): jxnEnds.add(p[0:i])
                for i in range(1,len(p)-2): jxnStarts.add(p[i:len(p)])



        obsSites = set(obsSites); jxnEnds=list(jxnEnds); jxnStarts=list(jxnStarts)

        jxnEnds.sort(key=lambda x: x[-1][1]); jxnStarts.sort()

        if len(jxnEnds) == 0 or len(jxnStarts) == 0: return

        jxnEnds = sorted([j for j in jxnEnds]); jxnStarts=sorted([j for j in jxnStarts])
        k=0
        for i in range(len(jxnEnds)):
            while k < len(jxnStarts) and jxnStarts[k][0][0] <= jxnEnds[i][-1][1]: k+=1
            for j in range(k,len(jxnStarts)):
                if j-k > self.maxJump: break
                if tuple((jxnEnds[i][-1][1],jxnStarts[j][0][0])) not in obsSites:
                    myCand =tuple([x for x in jxnEnds[i]]+[x for x in jxnStarts[j]])
                    if self.seqDist(myCand) > self.readlen:
                        novelSet.add(tuple([x for x in jxnEnds[i]]+[x for x in jxnStarts[j]]))
                        obsSites.add(tuple((jxnEnds[i][-1][1],jxnStarts[j][0][0])))
            
        self.novelJxns = sorted([s for s in novelSet])
            
            





    def findIntrons(self):
        introns=[]; maxlen = self.readlen -1 
        for i in range(len(self.exons)):
            if i > 0:
                EX=self.exons[i]; PREV=self.exons[i-1] 
                leftFlank  = max(PREV[1]-maxlen+1,self.start)
                rightFlank = min(EX[0] + maxlen-1, self.end)
                introns.append([(leftFlank,PREV[1]), (PREV[1]+1,EX[0]-1), (EX[0],rightFlank)])
        self.introns = introns



    def compilePairs(self):
        tmpPairs = [];  self.endTable = dd(list); self.startTable=dd(list)
        for p in self.knownJxns:
            tmpPairs += [(p[i-1][1],p[i][0]) for i in range(1,len(p))]
        self.validPairs = sorted([t for t in set(tmpPairs)])
        for t in self.validPairs:
            self.endTable[t[0]].append(t[1]); self.startTable[t[1]].append(t[0])





#######################################################################################################################################################################


                    

    def validOffsets(self):
        if self.name==None:
            return False
        else:
            self.checkOffsets()
            self.findMajorExons()
            self.findKnownJunctions()
            self.findNewJunctions()
            self.findIntrons()
            return True


    @staticmethod
    def uniquifyEnds(jxnList): 
        INTERNAL = set([]); MINSTARTS=dd(list); MAXENDS=dd(list); 
        for P in jxnList:
            INTERNAL.add(tuple([P[0][1]]+[tuple([t for t in P[1:len(P)-1]])]+[P[-1][0]]))
            for X in P:
                MINSTARTS[X[1]].append(X[0]); MAXENDS[X[0]].append(X[1])

        return list(set(sorted([ (((min(MINSTARTS[I[0]]),I[0]),) + I[1]  + ((I[2],max(MAXENDS[I[2]])),)) for I in INTERNAL])))



    @staticmethod
    def uniquifyMultis(jxnList):
        INTERNAL = set([]); MINSTARTS=dd(list); MAXENDS=dd(list); 
        for P in jxnList:
            
            MID_TUP=tuple([P[0][1]]+[tuple([t for t in P[1:len(P)-1]])]+[P[-1][0]])
            INTERNAL.add(MID_TUP)
            MINSTARTS[MID_TUP].append(P[0][0]); MAXENDS[MID_TUP].append(P[-1][-1])
            

        return list(set(sorted([ (((min(MINSTARTS[I]),I[0]),) + I[1]  + ((I[2],max(MAXENDS[I])),)) for I in INTERNAL])))
        
        


    @staticmethod
    def seqDist(jxnList):
        return sum([(x[1]-x[0])+1 for x in jxnList])





        
###########################################################################################################################################
