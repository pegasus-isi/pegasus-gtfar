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


def fileExtension(fname):
    return fname.split(".")[len(fname.split("."))-1]



def errorQuit(outString):
    sys.stderr.write(outString+'\n')
    sys.exit(2)


def listToString(myList,delims):
    


    if len(delims) == 3:
        
        topLevel = []
        for i in range(len(myList)):
            nextLevel = []
            for group in myList[i]:
                nextLevel.append(delims[0].join([str(item) for item in group]))
            topLevel.append(delims[1].join([element for element in nextLevel]))

        return delims[2].join([t for t in topLevel])

    print "BROKEN delims"
    sys.exit()
                

def locToCigar(spots):
    cigarStr = ''; spotList = [s for s in spots]; n=0

    while n+1<len(spotList):
        if spotList[n]+1 == spotList[n+1]:
            spotList.remove(spotList[n+1]); spotList.remove(spotList[n])
        else:
            n+=1
    for k in range(0,len(spotList),2):
        cigarStr+=str(spotList[k+1]-spotList[k]+1)+"M"
        if k+2< len(spotList):
            cigarStr+= str(spotList[k+2] - spotList[k+1] -1)+"N"
    return cigarStr
    

def cigarToLoc(loc,cigar):

    cigar= cigar.split("M")
    
    locs = [loc,loc + int(cigar[0]) -1]


    for c in cigar[1:len(cigar)-1]:
        c=c.split("N")
        locs.append(locs[-1] + int(c[0])+1)
        locs.append(locs[-1] + int(c[1])-1)
        
    return locs


def smartAvg(myList):
    if len(myList) == 0:
        return 0
    else:
        return sum(myList) / float(len(myList))

def smartDivide(num,denom):
    if denom ==0: return 0 
    else: return float(num)/denom




















