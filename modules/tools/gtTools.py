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

def seqCreate(seq,strand):
    newSeq=[]
    if strand == '-':
        baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        r=seq[-1::-1]
        for s in r:
            newSeq.append(baseComplement[s.capitalize()])
        return newSeq
    else:
        newSeq.extend([base.capitalize() for base in seq])
        return newSeq

def baseSwap(base):
    if base=="A": return 0
    if base=="C": return 1
    if base=="G": return 2
    if base=="T": return 3
    if base=="N": return 4
    
    if base==0: return 'A'
    if base==1: return 'C'
    if base==2: return 'G'
    if base==3: return 'T'
    if base==4: return 'N'

def fileExtension(fname):
    return fname.split(".")[len(fname.split("."))-1]



