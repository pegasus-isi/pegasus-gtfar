#!/usr/bin/env python


import sys
from MutationCand import *
from Sequence import *
from Utilities import *


from math import fabs


from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

class MutationFile:
    def __init__(self,mutationFile,minProb=0.6,minUrn=0.95):
        self.open = True
        self.minProb, self.minUrn = minProb, minUrn
        self.mFile = mutationFile
        self.fileHandle = open(mutationFile)
        self.nextLine()




###########################################################################################################################################################

    def nextLine(self):
        
        if self.open:
            
            self.line = self.fileHandle.readline().split()
            if len(self.line) < 1:
                self.open = False
                self.chr = None
                return 
            self.chr = self.line[0]
            self.pos = int(self.line[2])
            self.ref ,self.snp = self.line[3],self.line[4]
            self.cov, self.prob, self.urn= float(self.line[5]),float(self.line[6]), float(self.line[7])
###########################################################################################################################################################
