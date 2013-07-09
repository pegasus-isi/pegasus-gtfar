#!/usr/bin/env python


import sys
import difflib
import cPickle as pickle



#from gtTools import seqComplement


##########################################################################################################################################
##################################################  MULTI-GENE CLASS START  ##############################################################
##########################################################################################################################################

class GtfLine:
    def __init__(self,lp):
        line=lp.split();
        if len(line)==0:
            self.valid=False
            self.chr = 'NA'
        else:
            self.valid=True
            self.chr, self.db, self.type, self.start,self.end, self.strand= line[0], line[1],line[2],int(line[3]),int(line[4]),line[6]
            self.geneID, self.tranID, self.geneType, self.geneStatus, self.hugoName = line[9].split('"')[1], line[11].split('"')[1],line[13].split('"')[1], line[15].split('"')[1], line[17].split('"')[1]
