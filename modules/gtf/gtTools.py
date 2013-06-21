#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd


##########################################################################################################################################
##########################################################################################################################################
#####################################################  GT-FAR  FUNCTIONS #################################################################
##########################################################################################################################################
##########################################################################################################################################

#1) Complment of a sequence - input is either list OR string : returns string 

def seqComplement(seq):
    baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    p=""
    r=seq[-1::-1]
    for s in r:
        p+=baseComplement[s.capitalize()]
    return p

def fileExtension(fname):
    return fname.split(".")[len(fname.split("."))-1]


#class GtfLine:
#    def __init__(self,lp):
#        line=lp.split();
#        self.chr, self.db, self.type, self.start,self.end, self.strand= line[0], line[1],line[2],int(line[3]),int(line[4]),line[6]
#        self.geneID, self.tranID, self.geneType, self.geneStatus, self.hugoName = line[9].split('"')[1], line[11].split('"')[1],line[13].split('"')[1], line[15].split('"')[1], line[17].split('"')[1]

