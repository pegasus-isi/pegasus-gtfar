#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd
from tools.gtTools import seqComplement


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

class SamRead:
    def __init__(self,readList,keys,idx):
        
        self.ID         = readList[0][0]
        self.seq        = readList[0][9]
        self.qual       = readList[0][10]
        self.subs       = readList[0][11]
        self.strands    = [readList[i][1] for i in range(len(readList))]
        self.locations  = [readList[i][2] for i in range(len(readList))]
        self.geneIDs    = ["|".join([s for s in readList[i][2].split("|")[0:8]]) for i in range(len(readList))]
        self.positions  = [int(readList[i][3]) for i in range(len(readList))]
        self.mapTypes   = [readList[i][5] for i in range(len(readList))]
        self.locKeys    = keys
        self.index      = idx
        self.relocated  = False

    def relocate(self):
        self.geneLocs=[]; self.hgLocs=[]
        for i in range(len(self.locations)):
            DIST=len(self.seq)-1; POST=False
            for k in range(len(self.locKeys[i][0])):
                if self.positions[i]>=self.locKeys[i][0][k][0] and self.positions[i]<=self.locKeys[i][0][k][1]:
                    ### HERE ## 
                    tran_tuples, gene_tuples, hg_tuples = self.locKeys[i][0][k], self.locKeys[i][1][k], self.locKeys[i][2][k]

                    # DISCREPANCY COMES FROM A ZERO INDEXED GENOME/1 INDEXED GENES #
                    if tran_tuples[1] >= self.positions[i]+DIST:
                        posOffset=self.positions[i]-tran_tuples[0]
                        self.geneLocs.append([gene_tuples[0]+posOffset, gene_tuples[0]+posOffset+DIST])
                        if hg_tuples[0] < hg_tuples[1]:
                            self.hgLocs.append([hg_tuples[0]+posOffset-1, hg_tuples[0]+posOffset+DIST-1])
                        if hg_tuples[0] > hg_tuples[1]:
                            self.hgLocs.append([hg_tuples[0]-posOffset+1,hg_tuples[0]-posOffset-DIST+1])
                        break
                    else:
                        posOffset=self.positions[i]-tran_tuples[0]
                        self.geneLocs.append([gene_tuples[0]+posOffset,gene_tuples[1]+1])
                        if hg_tuples[0] < hg_tuples[1]:
                            self.hgLocs.append([hg_tuples[0]+posOffset-1,hg_tuples[1]])
                        else:
                            self.hgLocs.append([hg_tuples[0]-posOffset+1,hg_tuples[1]])
                        POST=True
                        DIST-=tran_tuples[1]-(tran_tuples[0]+posOffset-1)

                else:
                    if POST:
                        tran_tuples, gene_tuples, hg_tuples = self.locKeys[i][0][k], self.locKeys[i][1][k], self.locKeys[i][2][k]
                        if tran_tuples[1]>=tran_tuples[0]+DIST:
                            self.geneLocs[-1].extend((gene_tuples[0]+1,gene_tuples[0]+DIST))
                            if hg_tuples[0] < hg_tuples[1]:
                                self.hgLocs[-1].extend((hg_tuples[0],hg_tuples[0]+DIST-1))
                            else:
                                self.hgLocs[-1].extend((hg_tuples[0],hg_tuples[0]-DIST+1))
                            break
                        else:
                            self.geneLocs[-1].extend((gene_tuples[0]+1,gene_tuples[1]+1))
                            self.hgLocs[-1].extend((hg_tuples[0],hg_tuples[1]))
                            DIST-=(tran_tuples[1]-tran_tuples[0]+1)

        ## NOW TEST IF THE POSITIONS ARE UNIQ OR NOT ##

        self.hgUniq=True; self.multiAnnos=True; geneUniq=True; GLuniq=True
        hg1=self.hgLocs[0]; myChr=self.locations[0].split("|")[3]; g1=self.geneIDs[0]; myStrand=self.locations[0].split("|")[2]; samStrand=self.strands[0]; gL1=self.geneLocs[0];
                


        for i in range(len(self.locations)):
            if self.hgLocs[i]!=hg1 or self.locations[i].split("|")[3]!=myChr or self.locations[i].split("|")[2]!=myStrand or self.strands[i]!=samStrand:
                self.hgUniq=False
            if self.geneLocs[i]!=gL1:
                GLuniq=False
            if self.geneIDs[i]!=g1:
                geneUniq=False

        
        if self.hgUniq and geneUniq and GLuniq:
            
            self.geneFeatures=[]
            self.multiAnnos=False

            ## TOTALLY UNIQ ## 

            self.hgLoc = (samStrand,myChr,myStrand,tuple(hg1))
            self.geneLoc = (samStrand,g1,tuple(gL1))
            for i in range(len(self.locations)):
                self.geneFeatures.append(self.locations[i])

        elif not self.hgUniq:

            ## NON UNIQ ##

            self.multiLocs=[]; self.hgKey=dd(list); self.geneFeatureKey=dd(list)
            for i in range(len(self.locations)):
                tmpHgLoc = (self.strands[i],self.locations[i].split("|")[3],self.locations[i].split("|")[2],tuple(self.hgLocs[i]))
                if tmpHgLoc not in self.multiLocs:
                    self.multiLocs.append(tmpHgLoc)
                tmpGeneLoc=(self.strands[i],self.geneIDs[i],tuple(self.geneLocs[i]))
                if tmpGeneLoc not in self.hgKey[tmpHgLoc]:
                    self.hgKey[tmpHgLoc].append((self.strands[i],self.geneIDs[i],tuple(self.geneLocs[i])))
                self.geneFeatureKey[self.geneIDs[i]].append(self.locations[i])
        else:
            ## UNIQ BUT MULTIPLE ANNOTATIONS/ OR STRANDS IF NOT STRAND SPECIFIC ##
            self.geneFeatureKey=dd(list); self.annos = []
            self.hgLoc = (samStrand,myChr,myStrand,tuple(hg1))
            
            for i in range(len(self.locations)):
                tmpGeneLoc = (self.strands[i],self.geneIDs[i],tuple(self.geneLocs[i]))
                if tmpGeneLoc not in self.annos:
                    self.annos.append(tmpGeneLoc)
                self.geneFeatureKey[self.geneIDs[i]].append(self.locations[i])
        self.relocated=True

    def makeCigar(self,hgLocation):
        
        if self.hgUniq:
            scr='255'
        else:
            scr='1'

        if hgLocation[2]=="+":
            spots=hgLocation[3]
            samSeq=self.seq
            samStrand='0'
        else:
            spots=hgLocation[3][-1::-1]
            samSeq=seqComplement(self.seq)
            samStrand='16'

        startLoc=spots[0]; cigarStr=''
        
        for i in range(0,len(spots),2):
            cigarStr+=str(spots[i+1]-spots[i]+1)+"M"
            if i+2<len(spots):
                cigarStr+=str(spots[i+2]-spots[i+1]-1)+"D"
        
        return self.ID+'\t'+samStrand+'\t'+hgLocation[1]+'\t'+str(startLoc)+'\t'+scr+'\t'+cigarStr+'\t*\t0\t0\t'+samSeq+'\t'+'\t'+self.qual+'\t'+self.subs
        


    def visualizeBAM(self,samOut,SHOWAMBIG=False):
        if not self.relocated:
            print "Relocated Method Must be Called before Print Method"
            sys.exit()

        if self.index==1:
            self.writeSamHeader(samOut)

        if self.hgUniq:
            samOut.write('%s\n' % self.makeCigar(self.hgLoc))
        elif SHOWAMBIG:
            for m in self.multiLocs:
                samOut.write('%s\n' % self.makeCigar(m))


    def visualizeGeneLocation(self,readOut):
        
        if self.hgUniq:
            
            if not self.multiAnnos:
                ensgName=self.geneLoc[1].split("|")[0]
                spotStr="-".join([str(s) for s in self.geneLoc[2]])
                readOut.write("%s UNIQ SINGLE %s %s %s %s SEQ/QUAL: %s %s\n" % (self.ID,self.hgLoc[1],self.geneLoc[0],ensgName,spotStr,self.seq,self.qual))
            else:

                tmpAnnos=sorted([(anno[1].split("|")[0],anno[0],"-".join([str(s) for s in anno[2]])) for anno in self.annos])
                for t in tmpAnnos:
                    readOut.write("%s UNIQ MULTI %s %s %s %s SEQ/QUAL: %s %s\n" % (self.ID,self.hgLoc[1],t[1],t[0],t[2],self.seq,self.qual))
            
            
        else:
            for k in self.hgKey.keys():
                tmpAnnos=sorted([(anno[1].split("|")[0],anno[0],"-".join([str(s) for s in anno[2]])) for anno in self.hgKey[k]])
                for t in tmpAnnos:
                    if len(tmpAnnos)==1:
                        readOut.write("%s AMBIG SINGLE %s %s %s %s SEQ/QUAL: %s %s\n" % (self.ID,k[1],t[1],t[0],t[2],self.seq,self.qual))
                    else:
                        readOut.write("%s AMBIG MULTI %s %s %s %s SEQ/QUAL: %s %s\n" % (self.ID,k[1],t[1],t[0],t[2],self.seq,self.qual))
                
                    



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



    

