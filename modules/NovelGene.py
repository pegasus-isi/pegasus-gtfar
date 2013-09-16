#!/usr/bin/env python


import sys
from Sequence import *
from Utilities import *


from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class NovelGene:
    def __init__(self,mapFile,keyFile,prefix,SPECIES="HUMAN"):
        self.mapFile = mapFile; self.mapHandle = open(self.mapFile); self.open = True; self.offset = [0,0]; self.chr = None
        self.nextLine()
        self.chr = self.mapChr
        if not self.open:       errorQuit("EMPTY LOCATION FILE")
        
        

        if SPECIES == "HUMAN":
            self.chr2Num ={'chr1': 1, 'chr2': 2, 'chr3': 3,'chr4': 4, 'chr5': 5, 'chr6': 6,'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr10': 10, 'chr11': 11, 'chr12': 12,
                                'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18,'chr19': 19, 'chr20': 20, 'chr21': 21,'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrM': 25}
        elif SPECIES == "RHESUS":
            self.chr2Num = {'chr01': 1, 'chr02a': 2, 'chr02b': 3, 'chr03': 4, 'chr04': 5, 'chr05': 6,'chr06': 7, 'chr07': 8, 'chr08': 9,'chr09': 10, 'chr10': 11, 'chr11': 12, 'chr12': 13,
                    'chr13': 14, 'chr14': 15, 'chr15': 16, 'chr16': 17, 'chr17': 18, 'chr18':19 ,'chr19': 20, 'chrX': 21}
        self.keyHandle = open(keyFile)

        self.key = dd(list)
        
        for line in self.keyHandle:
            line = line.split()
            tmpdata = line[0].split("|"); tmppos=line[3].split("-")
            if tmpdata[-1][0:4]=="GENE":
                gene_id = "|".join(tmpdata[0:4])
                gene_chr = tmpdata[3]
                gene_tuple=int(tmppos[0]),int(tmppos[-1])
                self.key[gene_chr].append([gene_tuple,gene_id])


##########################################################################################################################################################

    def nextLine(self):
        if self.open:
            self.line = self.mapHandle.readline().split()
            if len(self.line)==11:
                self.mapChr,self.start,self.end,self.strand,spliceData = self.line[1],int(self.line[3]),int(self.line[4]),self.line[2],self.line[10].split("|")
               
                self.dist = self.end - self.start + 1
                if spliceData[0] == "INTERGENIC-SPLICED":
                    self.spliceSite = [int(x) for x in spliceData[1].split(":")[1].split("-")]
                else:
                    self.spliceSite = self.strand
            elif len(self.line)==0:
                self.mapGene = None
                self.mapChr  = None
                self.open = False
                #self.mapHandle.close()
            else:
                print "WTF"
###########################################################################################################################################################
           
    def search(self):
        self.geneTable = None
        self.hgStarts, self.hgSplice = [[0,0]],[None]
        ## ASSUMPTION -- APPROXIMATE READ LENGTH IS ABOUT 100 - so we want about 1/10 cov for 10x ##
        while self.chr == self.mapChr:

            if self.start - self.hgStarts[-1][0] < 900:
                self.hgStarts.append([self.start,self.end])
                self.hgSplice.append(self.spliceSite)
            else:
                self.hgStarts = [[self.start,self.end]]
                self.hgSplice = [self.spliceSite]
            
            if len(self.hgStarts) > 30:
                SPAN = self.hgStarts[-1][0] - self.hgStarts[0][0]
                self.rate = float(len(self.hgStarts)) / SPAN
                if self.rate > 0.01:
                    self.analyze()
                    self.hgStarts, self.hgSplice = [[0,0]],[None]
                else:
                    k=1; tmpRate = self.rate
                    while tmpRate < 0.01 and k+3<len(self.hgStarts):
                        k+=1
                        tmpRate = float(len(self.hgStarts[k::])) / (self.hgStarts[-1][0] - self.hgStarts[k][0])
                    self.hgStarts = self.hgStarts[k::]; self.hgSplice = self.hgSplice[k::]


            
            self.nextLine()
        self.chr = self.mapChr 



    def analyze(self):
        self.valid = False
        self.removeFromLeft()
        self.addFromRight()
        self.findBorders()
        if self.valid:
            self.genePrint()

    def genePrint(self):
        
        if self.geneTable == None:
            self.geneTable = sorted(self.key[self.chr])

        i =0 ; j =0
        while i < len(self.clusters): 
            j=0; print i 
            while j <(len(self.geneTable)):
                myCluster = self.clusters[i]
                if self.geneTable[j][0][1] < self.clusters[i][0]:
                    if j > 2: self.geneTable.remove(self.geneTable[j-1])
                    else: j+=1
                else:
                    if self.geneTable[j][0][0] > self.clusters[i][0]:
                        prevGene = self.geneTable[j-1]; postGene = self.geneTable[j]
                        prevDist = myCluster[0] - prevGene[0][1]
                        postDist = postGene[0][0] - myCluster[3]
                        if prevDist < postDist:
                            minDist = prevDist; geneName = prevGene[1]
                            if geneName.split("|")[2] == "-":
                                TYPE = "5-PRIME START"
                            else:
                                TYPE = "3-PRIME START"
                        else:
                            minDist = postDist;
                            geneName = postGene[1]
                            if geneName.split("|")[2] == "-":
                                TYPE = "3-PRIME END"
                            else:
                                TYPE = "5-PRIME START"
                        
                        print self.chr,myCluster[0],myCluster[3],myCluster[3]-myCluster[0],"edge/body/edge:",myCluster[0],myCluster[1],myCluster[2],myCluster[3],"supporting-reads:",myCluster[4][0],myCluster[4][1],myCluster[4][2],
                        print "+/-:",myCluster[5]["+"],myCluster[5]['-'],"Nearest-Neighbor:",minDist,geneName,TYPE
                        i+=1
                        if i == len(self.clusters): return 

                    else:
                        print "WEIRD CASE ---> WE ARE IN THE GENE"
                        print self.clusters[i],self.geneTable[j],self.geneTable[j-1],self.geneTable[j+1]
            
            if j == len(self.geneTable):
                break




    def removeFromLeft(self):
        for i in range(len(self.hgStarts)):
            newScr = float(len(self.hgStarts[i::])) / (self.hgStarts[-1][0] - self.hgStarts[i][0])
            if newScr *0.90 <= self.rate: 
                break 
            else:
                self.rate = newScr
        self.hgStarts = self.hgStarts[i::]; self.hgSplice = self.hgSplice[i::]
 

    def addFromRight(self):
 
        newScr  = self.rate
        while newScr >= (self.rate * 0.90) or newScr >= tmpRate or newScr>0.04: 
            tmpRate = newScr
            self.nextLine()
            if self.start - self.hgStarts[-1][0] < 1000:
                self.hgStarts.append([self.start,self.end])
                self.hgSplice.append(self.spliceSite)
            else:
                break
            newScr = float(len(self.hgStarts)) / (self.hgStarts[-1][0] - self.hgStarts[0][0])
        prevScr = float(len(self.hgStarts[0:len(self.hgStarts)-1])) / (self.hgStarts[-2][0] - self.hgStarts[0][0])
        if prevScr > newScr:
            self.hgStarts.pop(); self.hgSplice.pop()
            self.rate = prevScr
        else:
            self.rate = newScr

    
    def clusterContReads(self,contReads):
        contReads.sort()
        strand_dict = {"+": 0,"-":0}
        clusterStart = [[0,0,0,0,strand_dict],contReads[0][0]]; count = 0
        for c in contReads[1::]:
            if c[0][0] < contReads[-1][0][-1] + 30:
                count+=1
                strand_dict[c[1]]+=1
                clusterStart[-1][-1] = max(c[0][1],clusterStart[-1][-1])
            else:
                clusterStart[-1].extend([count,clusterStart[-1][0]-clusterStart[-1][0]+1,strand_dict])
                clusterStart.append([c][0])
                strand_dict = {"+": 0,"-":0}
                count = 0
        clusterStart[-1].extend([count,clusterStart[-1][0]-clusterStart[-1][0]+1,strand_dict])
       
        clusterStart.append([clusterStart[-1][1]+1000,clusterStart[-1][1]+1000,count,0,strand_dict])
        self.contClusters = clusterStart




    def findBorders(self):
        
        left = dd(int); right = dd(int); contReads = []; minC = None; maxC = None; spliceDict = dd(int)
        for i in range(len(self.hgSplice)):
            if self.hgSplice[i] == "+" or self.hgSplice[i] == "-":
                cont = self.hgStarts[i]
                if minC == None:
                    minC = cont[0]; maxC= cont[1]
                    if cont[0] < minC: minC = cont[0]
                    if cont[1] > maxC: maxC = cont[1]
                contReads.append([self.hgStarts[i],self.hgSplice[i]])
            else:
                spliceDict[tuple(self.hgSplice[i])]+=1

        if minC == None or maxC == None:
            return False
        self.clusterContReads(contReads)
        spliceList = sorted([ [spliceDict[k],k] for k in spliceDict.keys()],reverse=True) 
        FINAL_CLUSTERS = []
        for i in range(1,len(self.contClusters)-1):
            myCluster= self.contClusters[i]
            leftCands = []; rightCands = []
            for j in range(len(spliceList)):

                ################ ADD BOTH #########
                if spliceList[j][1][-1] < myCluster[0] and spliceList[j][0]>1:
                    tmpDist = myCluster[0] - spliceList[j][1][-1]
                    leftCands.append([tmpDist,spliceList[j][0],spliceList[j][1]])
                elif spliceList[j][1][0] > myCluster[1]:
                    tmpDist = spliceList[j][1][0] - myCluster[1]
                    rightCands.append([tmpDist,spliceList[j][0],spliceList[j][1]])
            if len(leftCands) > 0 and len(rightCands) > 0:
                L1 = sorted(leftCands)[0]
                R1 = sorted(rightCands)[0]
                FINAL_CLUSTERS.append([ L1[2][1],myCluster[0],myCluster[1],R1[2][0], (L1[1],myCluster[2],R1[1]),myCluster[4]])
            leftCands = []; rightCands = []

        if len(FINAL_CLUSTERS)>0:
            self.valid = True 
            self.clusters = FINAL_CLUSTERS
        




