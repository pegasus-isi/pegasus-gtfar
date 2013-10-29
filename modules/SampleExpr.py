#!/usr/bin/env python


import sys
from random import randrange
from random import choice
from collections import defaultdict as dd
from numpy import median
from numpy import mean

##########################################################################################################################################
#####################################################  GENE CLASS START   ################################################################
##########################################################################################################################################


class SampleExpr:
    def __init__(self,fname):
        self.name = fname
        self.raw  = dd(float)
        self.norm  = dd(float)
        self.rawCnts  = []



    def addGeneCnt(self,gene,cnt):

        self.raw[gene]=cnt
        self.rawCnts.append(cnt)

    def process(self,genes,cnt):
        myRatios = [];
        for g in genes:
            if self.raw[g] > 0: myRatios.append(self.raw[g] / (genes[g]**(1.0 / cnt)))
        
        self.medRatios  = median(myRatios)
        self.medianCnts = median(self.rawCnts)
        self.meanCnts   = mean(self.rawCnts)
        self.totals     = sum(self.rawCnts)

    def addResult(self,scr):
        self.phenotype = scr

    def normalize(self,method):
        if method == "MEDIAN":
            for g in self.raw.keys():
                self.norm[g] = self.raw[g] /float(self.medianCnts)
        if method == "MEAN":
            for g in self.raw.keys():
                self.norm[g] = self.raw[g] /float(self.meanCnts)
        if method == "NONE":
            for g in self.raw.keys():
                self.norm[g] = self.raw[g]
        if method == "RATIO":
            for g in self.raw.keys():
                self.norm[g] = self.raw[g] / float(self.medRatios) 

