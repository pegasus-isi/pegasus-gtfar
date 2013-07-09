#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd
from math import fabs
from ..tools.gtTools import *


##########################################################################################################################################
#####################################################  SAM READ  CLASS START ############################################################
##########################################################################################################################################

class SamRead:
    def __init__(self,readList,keys,SENSE):
        
        self.ID         = readList[0][0];  self.seq        = readList[0][9]
        self.qual       = readList[0][10]; self.subs       = readList[0][11]
        
        self.complement = None

        self.data    = [[s for s in readList[i][2].split("|")] for i in range(len(readList))]
        
        self.genomeStrand = [self.data[i][2] for i in range(len(self.data))]
        self.mapStrand    = [readList[i][1] for i in range(len(readList))]

        self.featurePos= [int(readList[i][3])-1 for i in range(len(readList))]

        self.sense = SENSE
        self.locKeys    = keys
        
        if keys != None:
            self.featureKeys=[keys[i][0]  for i in range(len(keys))]
            self.geneKeys=[keys[i][1]  for i in range(len(keys))]
            self.hgKeys=[keys[i][2]  for i in range(len(keys))]

        self.samStrings =[]; self.locStrings = []

##################################################################################################################################################################################
##################################################################################################################################################################################
#########################################################################  READ PROCESSES ########################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################



# 1)  seqComplement - Returns and creates the reverse complement of the read sequence 


    def seqComplement(self):
        if not self.complement:
            self.complement=""
            baseComplement={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
            for s in self.seq[-1::-1]:
                self.complement+= baseComplement[s.capitalize()]
        return self.complement


##################################################################################################################################################################################


# 2) Relocates the read - The mapping is initially in terms of the 'feature' (eg. Transcript, Intron, etc).  This method converts the feature position in gene and genome coordinates 


    def relocate(self):

        relocations=[]
        for i in range(len(self.featureKeys)):
            POST=False; DIST=len(self.seq)-1
            genePos=[]; hgPos =[]
            for j in range(len(self.featureKeys[i])):
                if POST or (self.featurePos[i] >= self.featureKeys[i][j][0] and self.featurePos[i] <= self.featureKeys[i][j][1]):
                    if not POST:
                        myOffset = self.featurePos[i] - self.featureKeys[i][j][0]

                    genePos.extend([self.geneKeys[i][j][0]+ myOffset, self.geneKeys[i][j][0] + myOffset + DIST  ])


                    if self.genomeStrand[i]=="+": hgPos.extend([ self.hgKeys[i][j][0] + myOffset, self.hgKeys[i][j][0] + myOffset + DIST ])
                    else:  hgPos.extend ( [ self.hgKeys[i][j][0] - myOffset, self.hgKeys[i][j][0] - myOffset - DIST ] )
                
                    if self.featureKeys[i][j][1] >= self.featurePos[i]+DIST and self.featureKeys[i][j][1] >= self.featureKeys[i][j][0]+DIST:
                        break
                    else:
                        genePos[-1]= self.geneKeys[i][j][1]; hgPos[-1]=self.hgKeys[i][j][1]
                        DIST -= (1 + self.featureKeys[i][j][1] - ( self.featureKeys[i][j][0] + myOffset ) )
                        POST=True; myOffset = 0
            relocations.append([(self.data[i][3],self.genomeStrand[i],hgPos),(self.data[i][0],self.mapStrand[i],genePos)])
        relocations.sort()
        self.locations = [relocations[0]]
        for r in relocations[1::]:
                
            if r[0]!=self.locations[-1][0]:        self.locations.append(r)
            elif r[1] not in self.locations[-1]:  self.locations[-1].append(r[1])



##################################################################################################################################################################################



# 3) Gathers the set of junctions which the alignemnt spans 



    def gatherJxns(self,minOverLap):

        self.junctionSpan = None
        if len(self.locations)>1 or len(self.locations[0][0][2])<3: return 
        else:
            genomeJumps=self.locations[0][0][2]; hJxns=[]; rJxns=[[self.locations[0][i][0]] for i in range(1,len(self.locations[0]))]
            
            for i in range(0,len(genomeJumps)-2,2):
                if fabs(genomeJumps[i+1]-genomeJumps[i]) > minOverLap and fabs(genomeJumps[i+3]-genomeJumps[i+2]) > minOverLap:
                        hJxns.append(str(genomeJumps[i+1])+'>'+str(genomeJumps[i+2]))
                        for j in range(len(rJxns)):
                            rJxns[j].append((str(self.locations[0][j+1][2][i+1])+'>'+str(self.locations[0][j+1][2][i+2])))

            if hJxns!=[]:
                mapType=set([]); k=0
                for d in self.data:
                    mapType.add(d[6][0:4])
                mapType=sorted(list(mapType))
                if   'TRAN' in mapType:     strStart='KNOWN_EXON_EXON_JXN|'
                elif 'ITRN' in mapType:     strStart='INRON_EXON_JXN|'
                elif mapType == ['GENE']:   strStart='NOVEL_EXON_EXON_JXN|'
                elif mapType == ['FLNK']:
                    if rJxns[0][1][0]<0:    strStart='5P_PREGENE_JUNCTION|'
                    else:                   strStart='3P_POSTGENE_JUNCTION|'
                
                strStart+=self.locations[0][0][0]+"|"+self.locations[0][0][1]+"|"+','.join(hJxns)+'|GENES|'
                for r in sorted(rJxns):
                    if k==0:
                        strStart+=r[0]+'|'+",".join([x for x in r[1::]])
                    elif k>2:
                        strStart+="&CONT..."
                        break
                    else:
                        strStart+='-&-'+r[0]+'|'+",".join([x for x in r[1::]])
                    k+=1;
                self.junctionSpan = strStart
        return
                        

##################################################################################################################################################################################



# 4) Converts location information into a str for the creation of a sam file (genome) and location file (genes)  


    def visualizeStrings(self):

        chrDict ={'chr1': 1, 'chr2': 2, 'chr3': 3,'chr4': 4, 'chr5': 5, 'chr6': 6,'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr10': 10, 'chr11': 11, 'chr12': 12, 
                'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18,'chr19': 19, 'chr20': 20, 'chr21': 21,'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrMT': 25}


    
        if len(self.locations)>1:
            samScr='1';   hgType="AMBIG"
        else:
            samScr='255'; hgType="UNIQ"

        for location in self.locations:
            genomeLoc = location[0];   gType="SINGLE";              geneLoc   = location[1::];    cigarStr=''
    
            if genomeLoc[1]=="+":
                samSeq=self.seq;       spots=genomeLoc[2];          samStrand='0'
            else:
                samSeq=self.seqComplement();  spots=genomeLoc[2][-1::-1];  samStrand='16'
                
            for i in range(0,len(spots),2):
                cigarStr+=str(spots[i+1]-spots[i]+1)+"M"
                if i+2< len(spots):
                    cigarStr+= str(spots[i+2] - spots[i+1] -1)+"D"

            self.samStrings.append(self.ID+'\t'+samStrand+'\t'+genomeLoc[0]+'\t'+str(spots[0])+'\t'+samScr+'\t'+cigarStr+'\t*\t0\t0\t'+samSeq+'\t'+'\t'+self.qual+'\t'+self.subs)

            for gene in geneLoc:
                spotStr=','.join([str(s) for s in gene[2]])

                if len(geneLoc)>1:
                    gType="MULTI"

                geneString=" ".join([self.ID,hgType,gType,genomeLoc[0],gene[1],gene[0],spotStr,"SEQ/QUAL",self.seq,self.qual,str(chrDict[genomeLoc[0]])])
                if geneString not in self.locStrings:
                    self.locStrings.append(geneString)
                
                 
#########################################################################################################################################################################################################


































































































                

                    



    def writeSamHeader(self,outfile):
        
        hgHeader=["@HD\tVN:0.1.5c\tSO:queryname",
                "@SQ\tSN:chr10\tLN:135534748",
                "@SQ\tSN:chr11\tLN:135006517",
                "@SQ\tSN:chr12\tLN:133851896",
                "@SQ\tSN:chr13\tLN:115169879",
                "@SQ\tSN:chr14\tLN:107349541",
                "@SQ\tSN:chr15\tLN:102531393",
                "@SQ\tSN:chr16\tLN:90354754",
                "@SQ\tSN:chr17\tLN:81195211",
                "@SQ\tSN:chr18\tLN:78077249",
                "@SQ\tSN:chr19\tLN:59128984",
                "@SQ\tSN:chr1\tLN:249250622",
                "@SQ\tSN:chr20\tLN:63025521",
                "@SQ\tSN:chr21\tLN:48129896",
                "@SQ\tSN:chr22\tLN:51304567",
                "@SQ\tSN:chr2\tLN:243199374",
                "@SQ\tSN:chr3\tLN:198022431",
                "@SQ\tSN:chr4\tLN:191154277",
                "@SQ\tSN:chr5\tLN:180915261",
                "@SQ\tSN:chr6\tLN:171115068",
                "@SQ\tSN:chr7\tLN:159138664",
                "@SQ\tSN:chr8\tLN:146364023",
                "@SQ\tSN:chr9\tLN:141213432",
                "@SQ\tSN:chr17_ctg5_hap1\tLN:1680829",
                "@SQ\tSN:chr4_ctg9_hap1\tLN:590427",
                "@SQ\tSN:chr6_apd_hap1\tLN:4622291",
                "@SQ\tSN:chr6_cox_hap2\tLN:4795372",
                "@SQ\tSN:chr6_dbb_hap3\tLN:4610397",
                "@SQ\tSN:chr6_mann_hap4\tLN:4683264",
                "@SQ\tSN:chr6_mcf_hap5\tLN:4833399",
                "@SQ\tSN:chr6_qbl_hap6\tLN:4611985",
                "@SQ\tSN:chr6_ssto_hap7\tLN:4928568",
                "@SQ\tSN:chrM\tLN:16572",
                "@SQ\tSN:chr11_gl000202_random\tLN:40104",
                "@SQ\tSN:chr17_gl000203_random\tLN:37499",
                "@SQ\tSN:chr17_gl000204_random\tLN:81311",
                "@SQ\tSN:chr17_gl000205_random\tLN:174589",
                "@SQ\tSN:chr17_gl000206_random\tLN:41002",
                "@SQ\tSN:chr18_gl000207_random\tLN:4263",
                "@SQ\tSN:chr19_gl000208_random\tLN:92690",
                "@SQ\tSN:chr19_gl000209_random\tLN:159170",
                "@SQ\tSN:chr1_gl000191_random\tLN:106434",
                "@SQ\tSN:chr1_gl000192_random\tLN:547497",
                "@SQ\tSN:chr21_gl000210_random\tLN:27683",
                "@SQ\tSN:chr4_gl000193_random\tLN:189790",
                "@SQ\tSN:chr4_gl000194_random\tLN:191470",
                "@SQ\tSN:chr7_gl000195_random\tLN:182897",
                "@SQ\tSN:chr8_gl000196_random\tLN:38915",
                "@SQ\tSN:chr8_gl000197_random\tLN:37176",
                "@SQ\tSN:chr9_gl000198_random\tLN:90086",
                "@SQ\tSN:chr9_gl000199_random\tLN:169875",
                "@SQ\tSN:chr9_gl000200_random\tLN:187036",
                "@SQ\tSN:chr9_gl000201_random\tLN:36149",
                "@SQ\tSN:chrUn_gl000211\tLN:166567",
                "@SQ\tSN:chrUn_gl000212\tLN:186859",
                "@SQ\tSN:chrUn_gl000213\tLN:164240",
                "@SQ\tSN:chrUn_gl000214\tLN:137719",
                "@SQ\tSN:chrUn_gl000215\tLN:172546",
                "@SQ\tSN:chrUn_gl000216\tLN:172295",
                "@SQ\tSN:chrUn_gl000217\tLN:172150",
                "@SQ\tSN:chrUn_gl000218\tLN:161148",
                "@SQ\tSN:chrUn_gl000219\tLN:179199",
                "@SQ\tSN:chrUn_gl000220\tLN:161803",
                "@SQ\tSN:chrUn_gl000221\tLN:155398",
                "@SQ\tSN:chrUn_gl000222\tLN:186862",
                "@SQ\tSN:chrUn_gl000223\tLN:180456",
                "@SQ\tSN:chrUn_gl000224\tLN:179694",
                "@SQ\tSN:chrUn_gl000225\tLN:211174",
                "@SQ\tSN:chrUn_gl000226\tLN:15009",
                "@SQ\tSN:chrUn_gl000227\tLN:128375",
                "@SQ\tSN:chrUn_gl000228\tLN:129121",
                "@SQ\tSN:chrUn_gl000229\tLN:19914",
                "@SQ\tSN:chrUn_gl000230\tLN:43692",
                "@SQ\tSN:chrUn_gl000231\tLN:27387",
                "@SQ\tSN:chrUn_gl000232\tLN:40653",
                "@SQ\tSN:chrUn_gl000233\tLN:45942",
                "@SQ\tSN:chrUn_gl000234\tLN:40532",
                "@SQ\tSN:chrUn_gl000235\tLN:34475",
                "@SQ\tSN:chrUn_gl000236\tLN:41935",
                "@SQ\tSN:chrUn_gl000237\tLN:45868",
                "@SQ\tSN:chrUn_gl000238\tLN:39940",
                "@SQ\tSN:chrUn_gl000239\tLN:33825",
                "@SQ\tSN:chrUn_gl000240\tLN:41934",
                "@SQ\tSN:chrUn_gl000241\tLN:42153",
                "@SQ\tSN:chrUn_gl000242\tLN:43524",
                "@SQ\tSN:chrUn_gl000243\tLN:43342",
                "@SQ\tSN:chrUn_gl000244\tLN:39930",
                "@SQ\tSN:chrUn_gl000245\tLN:36652",
                "@SQ\tSN:chrUn_gl000246\tLN:38155",
                "@SQ\tSN:chrUn_gl000247\tLN:36423",
                "@SQ\tSN:chrUn_gl000248\tLN:39787",
                "@SQ\tSN:chrUn_gl000249\tLN:38503",
                "@SQ\tSN:chrX\tLN:155270561",
                "@SQ\tSN:chrY\tLN:59373567"
                "@RG\tID:knowles.fastq\tSM:knowles.fastq",
                "@PG\tID:PerM\tVN:0.4.0"]

        for h in hgHeader: outfile.write("%s\n" % h)



    



##############################################################################################################






