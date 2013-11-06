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

##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class GeneExpr:
    def __init__(self,fname,data):
        self.name = fname
        self.data = data

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



        
