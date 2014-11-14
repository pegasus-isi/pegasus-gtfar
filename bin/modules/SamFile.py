#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd


#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class SamFile:
    def __init__(self,fileHandle,prefix): #,strictLevel):
        
        self.fileHandle = fileHandle
        self.prefix = prefix


        ### INFER FILE TYPE FROM FIRST LINE ### 
        self.myLine=self.fileHandle.readline().split()
        self.name,self.nextName = self.myLine[0],self.myLine[0]
        self.fileOpen= True  
        ### MAKE DICTIONARIES ###

        self.subcnts = dd(int)
        self.typeData = dd(int)
        self.senseData = dd(int)
        self.multiAnnos = dd(int)
        self.multiLocs = dd(int)
        self.geneData = dd(lambda: [0,0,0,0,])
        self.intronData = dd(lambda: [0,0,0,0,])
        self.classType = dd(lambda: [0,0,0,0])
        self.featureData = dd(lambda: [0,0,0,0])
        self.chrData = dd(lambda: [0,0,0,0])
        self.familyData = dd(lambda: [0,0,0,0])
        self.multiData = dd(lambda: [0,0])
        self.totalUniq = 0 
        self.multiGene = 0
        self.multiLoc  = 0 
        self.multiBoth = 0
        
    def getSamReadData(self):
        

        ### USE A MORE DATA CENTRIC WAY OF LOOKING AT THINGS ###

        readLines = [] 

        while self.name == self.nextName:
            readLines.append(self.myLine)
            try:
                self.myLine=self.fileHandle.readline().split()
                self.nextName = self.myLine[0] 
            except IndexError:
                self.fileOpen = False
                break
        
        if len(readLines[0][11].split(":i:")) > 1:
            subs,GT,TT,className,geneName,altName,geneFam,sense,readType = [x.split(":i:")[-1] for x in readLines[0][11::]]
            SEP=":i:"
        elif len(readLines[0][11].split(":Z:")) > 1:
            subs,GT,TT,className,geneName,altName,geneFam,sense,readType = [x.split(":Z:")[-1] for x in readLines[0][11::]]
            SEP=":Z:"
        else:
            print "ERROR"
            sys.exit()
        ### CLASS NAME IS THE INFORMATION ABOUT WHICH CLASS IT IS EG: ITRN:COORDS=152865354-152869580 ###
        ### GENE FAM IS INFO LIKE: protein_coding ###
        self.subcnts[subs]+=1
        self.typeData[readType]+=1
        if sense == "True":
            idx=0
            self.senseData[0]+=1
        else:
           idx=1
           self.senseData[1]+=1


        if len(readLines) == 1:
            self.totalUniq+=1
            if className.split(":")[0] in ["EXON","FILTER","KJXN","NJXN"]:
                self.geneData[geneName+','+altName][idx]+=1
            elif className.split(":")[0] in ["ITRN","IJXN"]:
                self.intronData[geneName+','+altName][idx]+=1
            #print className
            if className.split(":")[0] in ["EXON","KJXN","NJXN","ITRN"]:
                self.featureData[geneName+','+altName+"@"+className][idx]+=1

        #### FIX IT SO THE CLASS TYPE HAS THE GENE NAME WITH IT!!! ####
            self.classType[className.split(":")[0]][idx]+=1
            self.chrData[readLines[0][2]][idx]+=1
            self.familyData[geneFam][idx]+=1
        else:
            myGenes = set([x[15].split(SEP)[-1]+","+x[16].split(SEP)[-1] for x in readLines])
            myClasses = set([x[14].split(SEP)[-1] for x in readLines])
            myGeneClasses = set([x[15].split(SEP)[-1]+","+x[16].split(SEP)[-1]+"@"+x[14].split(SEP)[-1] for x in readLines])
            myFams = set([x[17].split(SEP)[-1] for x in readLines])
            myLocs = set([x[2]+","+x[3] for x in readLines])
            myChrs = set([x[2] for x in readLines])
            if len(myGenes) == 1:
                self.multiLoc+=1
                if className.split(":")[0] in ["EXON","FILTER","KJXN","NJXN"]:
                    self.geneData[geneName+','+altName][idx]+=1
                elif className.split(":")[0] in ["ITRN","IJXN"]:
                    self.intronData[geneName+','+altName][idx]+=1
                #self.geneData[geneName+','+altName][idx]+=1
                self.chrData[readLines[0][2]][idx]+=1
                
                self.familyData[geneFam][idx]+=1
                if len(myClasses) == 1:
                    if className.split(":")[0] in ["EXON","KJXN","NJXN","ITRN"]:
                        self.featureData[geneName+','+altName+"@"+className][idx]+=1
                    self.classType[className.split(":")[0]][idx]+=1
                else:
                    for c in myClasses:
                        if c.split(":")[0] in ["EXON","KJXN","NJXN","ITRN"]:
                            self.featureData[geneName+','+altName+"@"+c][idx+2]+=1.0/len(myClasses)
                        self.classType[c.split(":")[0]][idx+2]+=1.0/len(myClasses)
            elif len(myLocs) == 1:
                self.multiGene+=1
                
                if className.split(":")[0] in ["EXON","FILTER","KJXN","NJXN"]:
                    multiGene = "/".join([g for g in myGenes])
                    self.geneData[multiGene][idx]+=1
                self.chrData[readLines[0][2]][idx]+=1
                if len(myClasses) == 1:
                    if className.split(":")[0] in ["EXON","KJXN","NJXN","ITRN"]:
                        self.featureData[geneName+','+altName+"@"+className][idx]+=1
                    self.classType[className.split(":")[0]][idx]+=1
                else:
                    for c in myClasses:
                        self.classType[c.split(":")[0]][idx+2]+=1.0/len(myClasses)
                    for g in myGeneClasses:
                        if g.split(":")[0].split("@")[-1] in ["EXON","KJXN","NJXN","ITRN"]:
                            self.featureData[g][idx+2]+=1.0/len(myGenes)
                if len(myFams) == 1:
                    self.familyData[geneFam][idx]+=1
                else:
                    for f in myFams:
                        self.familyData[f][idx]+=1.0/len(myFams)
                for g in myGenes:
                    self.geneData[g][idx+2]+=1.0/len(myGenes)
        #            for j in myGenes:
         #               self.geneData[g][4].add(j)
                tmpMultiGene = "/".join(sorted([g for g in myGenes]))
                self.multiData[tmpMultiGene][0]+=1
            else:
                self.multiBoth+=1
                for c in myChrs:
                    self.chrData[c][idx+2]+=1.0/len(myChrs)
                for g in myGeneClasses:
                    if g.split(":")[0].split("@")[-1] in ["EXON","KJXN","NJXN","ITRN"]:
                        self.featureData[g][idx+2]+=1.0/len(myClasses)
                for c in myClasses:
                    #self.featureData[geneName+','+altName+"@"+c][idx+2]+=1.0/len(myClasses)
                    self.classType[c.split(":")[0]][idx+2]+=1.0/len(myClasses)
                for f in myFams:
                    self.familyData[f][idx+2]+=1.0/len(myFams)
                for g in myGenes:
                    if className.split(":")[0] in ["EXON","FILTER","KJXN","NJXN"]:
                        self.geneData[g][idx+2]+=1.0/len(myGenes)
                tmpMultiGene = "/".join(sorted([g for g in myGenes]))
                self.multiData[tmpMultiGene][0]+=1
                
        self.name = self.nextName



    #self.writeTable = {self.readLen: open(self.prefix+'_full.fastq',"w"), 0: open(self.prefix+'_reject.fastq',"w") }
    def printResults(self):
        geneOut = open(self.prefix+'.gene.cnts','w')
        multiGeneOut = open(self.prefix+'.overlapGene.cnts','w')
        for g in self.geneData:
            tmpData=self.geneData[g]
            #if len(tmpData[4]) ==0: multiAnnos = "None"
            #else:                   multiAnnos = "/".join([s for s in tmpData[4]])
            #if len(tmpData[5]) ==0: multiLocs = "None"
            #else:                   multiLocs = "/".join([s for s in tmpData[5]])

            if len(g.split("/"))==1:
                if tmpData[0]+tmpData[1]>0:
                    geneOut.write('%s %s %s\n' % (g,tmpData[0],tmpData[1]))
            else:
                if tmpData[0]+tmpData[1]>0:
                    multiGeneOut.write('%s %s %s\n' % (g,tmpData[0],tmpData[1]))

        geneOut.close()
        multiGeneOut.close()

        featOut = open(self.prefix+'.feature.cnts','w')
        for c in self.featureData:
            featOut.write('%s %s %s\n' % (c,self.featureData[c][0],self.featureData[c][1])) 
            #" ".join([str(s) for s in self.featureData[c]])))
        featOut.close()
        
        multiOut = open(self.prefix+'.ambiguousGenes.cnts','w')
        
        for c in self.multiData:
            multiOut.write('%s %s\n' % (c," ".join([str(s) for s in self.multiData[c]])))





        summaryOut = open(self.prefix+'.summary.out','w')
        #summaryOut = sys.stdout

        summaryOut.write('Total-Reads: %s Unique %s MultiAnnotated %s RepetitiveGene %s Ambiguous %s\n' % (self.totalUniq+self.multiGene+self.multiLoc+self.multiBoth,self.totalUniq,self.multiGene,self.multiLoc,self.multiBoth ))
        summaryOut.write('Sense/AntiSense: %s %s\n' % (self.senseData[0],self.senseData[1]))
        
        
        subs_seen = sorted([int(s) for s in self.subcnts.keys()])
        classes_seen = sorted([c for c in self.classType.keys()])
        familys_seen = sorted([f for f in self.familyData.keys()])
        types_seen = sorted([f for f in self.typeData.keys()])
        chrs_seen = sorted([c for c in self.chrData.keys()])
        
        
        summaryOut.write('Substitutions')
        for s in subs_seen: summaryOut.write(' %s %s' % (s,self.subcnts[str(s)]))
        summaryOut.write('\n')

        summaryOut.write('FeatureTypes')
        for c in classes_seen:  summaryOut.write(' %s %s' % (c," ".join([str(s) for s in self.classType[c]])))
        summaryOut.write('\n')

        summaryOut.write('FamilyTypes')
        for c in familys_seen:  summaryOut.write(' %s %s' % (c," ".join([str(s) for s in self.familyData[c]])))
        summaryOut.write('\n')
        
        summaryOut.write('readTypes')
        for c in types_seen:  summaryOut.write(' %s %s' % (c,self.typeData[c]))
        summaryOut.write('\n')

        summaryOut.write('Chromosomes\n')
        for c in chrs_seen:  summaryOut.write('%s %s\n' % (c," ".join([str(s) for s in self.chrData[c]])))
            
            

            
    def dataFailure(self,msg="FOO"):
        sys.stderr.write(msg+"\n")
        sys.exit()



            
