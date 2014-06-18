#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd
from scipy import stats

#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class CntFiles:
    def __init__(self,files,keys,ktype,prefix): #,strictLevel):


        self.METHOD="UNIQUE_SENSE"
        if ktype == 'binary':
            self.groupDict = {}
            self.flist = []
            for line in open(keys):
                line=line.split()
                if line[1] == '0' or line[1] == '1':
                    self.groupDict[line[0]]=int(line[1])
                    self.flist.append(line[0])
                    self.keyType = 'binary'
                else:
                    sys.stderr.write("Warning: Non-binary key supplied - assuming continuous variable\n")
                    self.keyType = 'continuous'
                    self.groupDict[line[0]] = float(line[1])
                    self.flist.append(line[0])


        self.cnt_dict = dd(lambda: dd(lambda: [0,0,0,0]))
        self.total_cnts = dd(float)
        for f in files:
            if f not in self.flist:
                print "ERROR XYZ"
                sys.exit()
            for line in open(f):
                line = line.split()
                self.cnt_dict[line[0]][f] = [float(line[1]),float(line[2]),float(line[3]),float(line[4])]
                self.total_cnts[f] = sum(self.cnt_dict[line[0]][f])
                

        self.prefix = prefix


    def analyze(self):
       

        METHOD = "UNIQUE_SENSE"
        self.gene_stats = []
        if self.keyType == 'binary':
            for gene in self.cnt_dict.keys():
                INIT=True
                CNT_DATA = [[],[]]
                for sample in self.flist:
                    if INIT:
                        first_sample = self.total_cnts[sample]
                        INIT=False

                    if sample in self.cnt_dict[gene].keys():
                        tmp_cnts = [(first_sample* x) / (self.total_cnts[sample]) for x in  self.cnt_dict[gene][sample]]
                    else:
                        tmp_cnts = [0.0,0.0,0.0,0.0]

                    if self.METHOD == "UNIQUE_SENSE":
                        CNT_DATA[self.groupDict[sample]].append(tmp_cnts[0])
                
                mean0 = sum(CNT_DATA[0]) / len(CNT_DATA[0])
                mean1 = sum(CNT_DATA[1]) / len(CNT_DATA[1])
                if mean0 > mean1: FC = mean0/(mean1+0.01)
                else:             FC = -1*(mean1/(mean0+0.01))

                tval = stats.ttest_ind(CNT_DATA[0],CNT_DATA[1])[1]
            
                self.gene_stats.append([tval,FC,mean0,mean1,gene,CNT_DATA[0],CNT_DATA[1]])
                #print self.gene_stats
        elif self.keyType == 'continuous':
            
            for gene in self.cnt_dict.keys():
                INIT=True
                phenos = []; values = []
                for sample in self.flist:
                    if INIT:
                        first_sample = self.total_cnts[sample]
                        INIT=False
                    phenos.append(self.groupDict[sample])
                    if sample in self.cnt_dict[gene].keys():
                        tmp_cnts = [(first_sample* x) / (self.total_cnts[sample]) for x in  self.cnt_dict[gene][sample]]
                    else:
                        tmp_cnts = [0.0,0.0,0.0,0.0]
                    
                    if self.METHOD == "UNIQUE_SENSE":
                        values.append(tmp_cnts[0])

                R = stats.pearsonr(values,phenos)[0]
                S = stats.spearmanr(values,phenos)[0]
                    
                tmp_data = sorted([(values[i],phenos[i]) for i in range(len(values))])

                self.gene_stats.append([R*R,R,S,gene,tmp_data])
                    
                

    def printResults(self):
        if self.keyType == 'binary':
            for g in sorted(self.gene_stats):
                print g[0],g[1],g[2],g[3],g[4],"|"," ".join([str(s) for s in g[5]]),"|"," ".join([str(s) for s in g[6]])
        if self.keyType == 'continuous':
            self.gene_stats.sort(reverse=True)
            for g in self.gene_stats:
                print g[0],g[1],g[2],g[3]
                 

            

            
    def dataFailure(self,msg="FOO"):
        sys.stderr.write(msg+"\n")
        sys.exit()



            
