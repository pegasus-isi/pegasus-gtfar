#!/usr/bin/env python


import sys
import difflib



##########################################################################################################################################
##################################################  MULTI-GENE CLASS START  ##############################################################
##########################################################################################################################################

class ProgressBar:
    def __init__(self,fHandle,initString,dot,finish,hitCnt,VERBOSE):

        self.fHandle = fHandle; self.dot = dot; self.finish = finish; self.hitCnt = hitCnt; self.VERBOSE=VERBOSE
        if self.VERBOSE: self.fHandle.write("\n"+initString)
        self.cnt = 0 
        
    def increment(self):
        self.cnt += 1
        if self.cnt % self.hitCnt == 0 and self.VERBOSE: self.fHandle.write(self.dot)


    def complete(self):
        if self.VERBOSE: self.fHandle.write(self.finish+"\n")
