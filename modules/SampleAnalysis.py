#!/usr/bin/env python


import sys
from random import randrange
from random import choice
from collections import defaultdict as dd
from numpy import median
from numpy import mean
from scipy import stats
from math import fabs
from math import isnan
import itertools
from math import log
##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class SampleData:
    def __init__(self):
        self.members = []
    
    def addMember(self,member):
        self.members.append(member)

    def mergeGenes(self):
        self.genes = list(set(itertools.chain(*[x.geneList for x in self.members])))

    def pairwiseCorrelation(self,idx1,idx2):
        List1, List2 = [], [] 
        for gene in self.genes:
            A=self.members[idx1].cnts.exonic(gene,"TOTAL"); B=self.members[idx2].cnts.exonic(gene,"TOTAL")
            if A>0 and B>0:
                #List1.append(log(A+.00001))
                #List2.append(log(B+.00001))
                List1.append(log(A))
                List2.append(log(B))
                #List1.append(self.members[idx1].cnts.exonic(gene,"UNIQUE"))
                #List2.append(self.members[idx2].cnts.exonic(gene,"UNIQUE"))
            #print gene,self.members[idx1].cnts.exonic(gene,"UNIQUE"),self.members[idx2].cnts.exonic(gene,"UNIQUE")
        CORR = stats.pearsonr(List1,List2)[0]
        print CORR
       



class GeneCnts:
    def __init__(self):
        self.rawCnts = dd(lambda: [0,0,0,0,0,0])
        self.totals  = [0,0,0,0,0,0]
        
    def add(self,x):
        myName=x[0]
        myCnts=[int(y) for y in [x[3],x[4],x[6],x[7],x[9],x[10]]]
        self.rawCnts[myName] = [myCnts[i]+self.rawCnts[myName][i] for i in range(len(myCnts))]
        self.totals          = [myCnts[i]+self.totals[i] for i in range(len(myCnts))]
       
    def exonic(self,gene,TYPE="UNIQUE",NORMALIZE=None):
       
        if not NORMALIZE:
            if TYPE=="UNIQUE":
                return self.rawCnts[gene][0]
            elif TYPE=="AMBIG":
                return self.rawCnts[gene][1]
            elif TYPE=="TOTAL":
                return self.rawCnts[gene][0]+self.rawCnts[gene][1]


        
        


    

         


class Sample:
    def __init__(self,fname,INSTRUCTION):
        #self.name = fname
        #self.cntData = dd(lambda: [0,0,0,0,0,0])
        if INSTRUCTION=="LOAD_GENE_CNTS":
            self.name=fname.split("_gene.cnts")[0]
            self.cnts = GeneCnts()
            self.geneList = []
            for line in open(fname):
                line=line.split()
                self.geneList.append(line[0])
                self.cnts.add(line) 
            


    def calculate_continuous_correlation(self):


        phenotypes  = []
        genotypes   = []

        self.data.sort()

        for x in self.data:
            try:
                
                phenotypes.append(float(x[0]))
                genotypes.append(float(x[1]))
            except ValueError:
                    print "DOES NOT LOOK CONT"
                    sys.exit()
        if len([g for g in genotypes if g>0])>len(phenotypes)/2.0: 
            CORR = stats.pearsonr(genotypes,phenotypes)[0]
            print self.name,CORR,fabs(CORR)
        else:
            print self.name, 0 ,0 



        
