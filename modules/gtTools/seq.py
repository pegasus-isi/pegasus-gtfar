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


def revComp(seq):
    baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    if type(seq) == str:
        newOut=''
        for s in seq[-1::-1]:
            newOut+=baseComplement[s.capitalize()]
    elif type(seq) == list:
        newOut = []
        for s in seq[-1::-1]:
            newOut.append(baseComplement[s.capitalize()])
    else:
        print 'WRONG TYPE'
        sys.exit()
    return newOut



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

def trueStrand(s1,s2,spots):

    if s1 == s2:
        strand = '+'
        if spots[0] > spots[1]:
            spots = spots[-1::-1]
    else:
        strand = '-'
        if spots[0] < spots[1]:
            spots = spots[-1::-1]

    return strand,spots


def sam2Map(item):
    if type(item) == str:
        if item == '16': return '-'
        if item == '0':  return '+'


def cigar2List(loc,cigar):
    
    cigar=cigar[0:len(cigar)-1]
    cigar=cigar.split("M")
    init=int(cigar[0])
    mapList=[loc,init+loc-1]
    for pair in cigar[1::]:
        pair=pair.split("N")
        newTuple=[mapList[-1]+int(pair[0])+1,mapList[-1]+int(pair[0])+int(pair[1])]
        mapList=mapList+newTuple
    return mapList 



def diversity(baseList):
    total = sum(baseList)
    return sum(baseList) - total / (total)



def fileExtension(fname):
    return fname.split(".")[len(fname.split("."))-1]



