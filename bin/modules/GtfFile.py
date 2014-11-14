#!/usr/bin/env python

import os

from ToolSet import errorQuit, listToString
from GtfLine import *
from GtfGene import *
from GtfFilters import *
from collections import defaultdict as dd

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################

class GtfFile:
    def __init__(self, fileHandle, prefix=None, readlen=None, filterType=None, findCands=True):
        try:    self.fName = open(fileHandle)
        except TypeError:   errorQuit("A GTF FILE IS NOT SUPPLIED")

        self.prefix, self.readLen,self.findCands = prefix, readlen, findCands

        self.genes, self.seq = [], []
        self.open = True
        self.minLen, self.maxLen = 35, 400
        self.SNPCANDS = False

        if self.prefix != None:
            dir_name = os.path.dirname(self.prefix)
            f_name = os.path.basename(self.prefix)

            if dir_name != '' and not os.path.exists(dir_name):
                os.makedirs(dir_name)

            features = self.prefix + '_features.fa' if f_name else self.prefix + 'FEATURES.fa'
            chrs = self.prefix + '_chrs.fa' if f_name else self.prefix + 'GENOME.fa'
            jxnCands = self.prefix + '_jxnCands.fa' if f_name else self.prefix + 'SPLICES.fa'
            gene = self.prefix + '_geneSeqs.fa' if f_name else self.prefix + 'GENE.fa'

            self.featureFile = open(features, 'w')
            self.chrFile = open(chrs, 'w')
            self.headerFile = open(self.prefix + "_headers.txt", "w")
            self.headerFile.write("%s\n" % "\t".join(["@HD", "VN:0.1.5c", "SO:queryname"]))
            self.headerFile.write("%s\n" % "\t".join(["@SQ", "SN:chrR", "LN:10000"]))

            if self.findCands:
                self.geneFile = open(gene, 'w')
                self.candFile = open(jxnCands, 'w')

        if filterType == "HUMAN":
            for name, seq in GtfFilters(filterType).seqs: self.featureFile.write(
                "%s:0-%s\n%s\n" % (name, len(seq), seq))

        tmpLine = self.fName.readline().strip()
        while tmpLine[0] == "#":
            tmpLine = self.fName.readline().strip()
        self.line = GtfLine(tmpLine)
        self.chr = self.line.chr



############################################################################################################################################
######################################## Annotation Methods ################################################################################
############################################################################################################################################


    def loadGenesOnChromosome(self):
        while self.chr == self.line.chr:
            self.gene = GtfGene(self.line,self.readLen)
            while self.line.valid and self.line.geneID ==  self.gene.name:
                self.gene.addGtfLine(self.line)
                self.line = GtfLine(self.fName.readline())
            if self.gene.validOffsets():
                self.genes.append(self.gene)


    def startNextChromosome(self):
        if self.line.chr == 'NA':
            self.chr =  'NA'
            self.open = False
        else:
            self.chr = self.line.chr; self.genes=[]; self.seq=[]; self.geneKey={}


    def addFasta(self,filePath):
        try:
            c=open(filePath)
        except IOError:
            if len(filePath.split("/")) > 1:
                failChr = filePath.split("/")[-1]
            else: failChr = filePath
            errorQuit("Error: "+failChr+" is referenced by the suppled gtf-file but not found in the supplied chromosome path: "+filePath)
        fChr = c.readline().strip().split()[0].split(">")[1]
        #if self.chr != c.readline().strip().split(">")[1]:
        if self.chr != fChr:
            print self.chr,fChr

            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])


    def printGenesOnChromosome(self,TYPE='ALL'):
        geneTuples=[]
        for gene in self.genes:
            geneTuples.append((gene.start-1,gene.end))
            geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            geneInfo = gene.name+"|"+gene.hugo+"|"+gene.type+"|"+gene.chr+"|"+gene.strand+"|"
            if gene.end - gene.start > self.minLen:
                if self.findCands:
                    self.geneFile.write(">%s\n%s\n" % ( geneInfo+"GENE:"+str(gene.start)+"-"+str(gene.end), geneSeq))
                    for N in gene.novelJxns: self.candFile.write(">%s\n%s\n" % (geneInfo+"NJXN:"+listToString(N,["|","-"]),"".join([geneSeq[n[0]-gene.start:(n[1]-gene.start)+1] for n in N])))

            for start,end in gene.exons:    self.featureFile.write(">%s\n%s\n" % (geneInfo+"EXON:"+str(start)+"-"+str(end),geneSeq[start-gene.start:(end-gene.start)+1]))
            for I in gene.introns:   self.featureFile.write(">%s\n%s\n" % (geneInfo+"INTRON:"+listToString(I,["|","-"]),"".join([geneSeq[i[0]-gene.start:(i[1]-gene.start)+1] for i in I])))
            for J in gene.knownJxns: self.featureFile.write(">%s\n%s\n" % (geneInfo+"KJXN:"+listToString(J,["|","-"]),"".join([geneSeq[j[0]-gene.start:(j[1]-gene.start)+1] for j in J])))


        for g in geneTuples:
            if (g[1]-g[0]) > (self.readLen+1)*2:
                self.seq[g[0]+self.readLen:g[1]-self.readLen]=["N" for i in range((g[1]-self.readLen)-(g[0]+self.readLen))]

        chrLen=len(self.seq)


        self.headerFile.write("%s\n" % "\t".join(["@SQ","SN:"+self.chr,"LN:"+str(chrLen)]))
        k=0;BUFFER=200
        self.chrFile.write("%s\n" % (">"+self.chr))
        while True:
            self.chrFile.write("%s\n" % "".join([b for b in self.seq[BUFFER*k:BUFFER*(k+1)]]))
            k+=1
            if k*BUFFER > chrLen: break








############################################################################################################################################
############################################################################################################################################
############################################################################################################################################


############################################################################################################################################
####################################################### CANDIDATE CODE  ####################################################################
############################################################################################################################################

    def addCandidates(self,candidate_file,cand_type):
        if cand_type=="SNPS":
            self.SNPCANDS = True
            self.snpcands = dd(lambda: dd(list))
            self.candFile    =   open(self.prefix+'_rddCands.fa','w')
            self.referenceFile    =   open(self.prefix+'_rddReference.txt','w')
            for line in open(candidate_file):
                line=line.split()
                self.snpcands[line[0]][int(line[1])]=[line[2::]]


    def printSnpCandsOnChromosome(self,TYPE='ALL'):
        geneTuples=[]
        if self.SNPCANDS and self.snpcands[self.chr] != []:
            snp_key = self.snpcands[self.chr]
            snp_locs = sorted(snp_key.keys())
            k=0
            if len(snp_locs) == 0: return

            for gene in self.genes:
                geneTuples.append((gene.start-1,gene.end))
                geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
                geneInfo = gene.name+"|"+gene.hugo+"|"+gene.type+"|"+gene.chr+"|"+gene.strand+"|"

                while snp_locs[k] < gene.start:
                    myLoc = snp_locs[k]
                    myKey = snp_key[myLoc]
                    refAllele = myKey[0][0]
                    altAlleles = myKey[0][1::]

                    intergenic_seq = ["".join(self.seq[(myLoc-1)-self.readLen:myLoc-1])]+[self.seq[myLoc-1]]+["".join(self.seq[myLoc:myLoc+self.readLen])]
                    intergenic_loc = [["intergenic"], myLoc,"NULL",[myLoc-self.readLen,myLoc+self.readLen+1],intergenic_seq]
                    self.referenceFile.write("%s %s %s %s NA NA NA " % (self.chr,myLoc,refAllele,"/".join(altAlleles)))
                    self.referenceFile.write("INTERGENIC NA NA %s\n" % "-".join([str(s) for s in intergenic_loc[3]]))
                    self.candFile.write(">%s_%s_NA_NA_NA_%s_%s_" % (self.chr,myLoc,"intergenic","-".join([str(s) for s in intergenic_loc[3]])))
                    self.candFile.write("REF-ALLELE_%s\n%s\n" %(refAllele,"".join([intergenic_seq[0],refAllele,intergenic_seq[-1]])))
                    n=0
                    for a in altAlleles:
                        n+=1
                        self.candFile.write(">%s_%s_NA_NA_NA_%s_%s_" % (self.chr,myLoc,"intergenic","-".join([str(s) for s in intergenic_loc[3]])))
                        self.candFile.write("ALT-ALLELE-%s_%s\n%s\n" %(str(n),a,"".join([intergenic_seq[0],a,intergenic_seq[-1]])))
                    k+=1
                    if k == len(snp_locs): return

                while snp_locs[k] <= gene.end:
                    myLoc = snp_locs[k]
                    myKey = snp_key[myLoc]
                    refAllele = myKey[0][0]
                    altAlleles = myKey[0][1::]
                    exonic_locs = []
                    prev_seqs = []

                    for t in gene.transcripts:
                        tran_exons = sorted(t[1])
                        for i in range(len(tran_exons)):
                            if myLoc >= tran_exons[i][0] and myLoc <= tran_exons[i][1]:
                                if myLoc - tran_exons[i][0] > self.readLen:
                                    left_side = [(myLoc-self.readLen,myLoc-1)]
                                else:
                                    left_side = [(tran_exons[i][0],myLoc-1)]
                                    remainder = self.readLen - (((myLoc-1) - tran_exons[i][0]) + 1)
                                    j=i-1
                                    while j >= 0 and (tran_exons[j][1] - tran_exons[j][0]) <= remainder:
                                        remainder -= ((tran_exons[j][1] - tran_exons[j][0]) + 1 )
                                        left_side.append(tran_exons[j])
                                        j-=1
                                    left_side.append(((tran_exons[j][1]-remainder)+1,tran_exons[j][1]))
                                    left_side.reverse()
                                if tran_exons[i][1] - myLoc > self.readLen:
                                    right_side = [(myLoc+1,myLoc+self.readLen)]
                                else:
                                    remainder = self.readLen - ((tran_exons[i][1] - (myLoc+1)) + 1)
                                    right_side = [(myLoc+1,tran_exons[i][1])]
                                    j=i+1
                                    while j < len(tran_exons) and tran_exons[j][1] - tran_exons[j][0] <= remainder:
                                        remainder -= ((tran_exons[j][1] - tran_exons[j][0]) + 1 )
                                        right_side.append(tran_exons[j])
                                        j+=1
                                    if j<len(tran_exons): right_side.append((tran_exons[j][0],(tran_exons[j][0] + remainder) - 1))

                                if self.seq[myLoc-1] != myKey[0][0]:
                                    print "WTF",self.seq[myLoc],myKey[0][0],self.seq[myLoc:myLoc+5]
                                if left_side[-1][1] +1 == myLoc and myLoc == right_side[0][0] -1:
                                    minDist = min(myLoc - left_side[-1][0],right_side[0][1]-myLoc)
                                    exonic_seq = ["".join(["".join(self.seq[x[0]-1:x[1]]) for x in left_side])]+[self.seq[myLoc-1]]+["".join(["".join(self.seq[x[0]-1:x[1]]) for x in right_side])]
                                    if exonic_seq not in prev_seqs:
                                        prev_seqs.append(exonic_seq)
                                        left_side[-1]=(left_side[-1][0],right_side[0][1])
                                        exonic_info = [[t[0]],myLoc,minDist,left_side+right_side[1::],exonic_seq]
                                        exonic_locs.append(exonic_info)
                                    else:
                                        for info in exonic_locs:
                                            if exonic_seq == info[-1]:
                                                info[0].append(t[0])
                                else:
                                    print "FIX"
                                    sys.exit()

                                break
                    intronic_seq = ["".join(self.seq[(myLoc-1)-self.readLen:myLoc-1])]+[self.seq[myLoc-1]]+["".join(self.seq[myLoc:myLoc+self.readLen])]
                    intronic_loc = [["intronic"], myLoc,"NULL",[myLoc-self.readLen,myLoc+self.readLen+1],intronic_seq]
                    if exonic_locs != []:
                        m=0
                        transcript_groupings = []
                        jxn_groups = []
                        for e in exonic_locs:
                            m+=1
                            transcript_groupings.append(",".join(e[0]))
                            jxn_groups.append("&".join([str(x[0])+"-"+str(x[1]) for x in e[3]]))
                            self.candFile.write(">%s_%s_%s_%s_%s_%s_%s_" % (self.chr,myLoc,gene.name,gene.hugo,gene.strand,"exonic"+str(m),jxn_groups[-1]))
                            self.candFile.write("REF-ALLELE_%s\n%s\n" %(refAllele,"".join([e[-1][0],refAllele,e[-1][-1]])))
                            n=0
                            for a in altAlleles:
                                n+=1
                                self.candFile.write(">%s_%s_%s_%s_%s_%s_%s_" % (self.chr,myLoc,gene.name,gene.hugo,gene.strand,"exonic"+str(m),jxn_groups[-1]))
                                self.candFile.write("ALT-ALLELE-%s_%s\n%s\n" %(str(n),a,"".join([e[-1][0],a,e[-1][-1]])))
                        self.referenceFile.write("%s %s %s %s %s %s %s " % (self.chr,myLoc,refAllele,"/".join(altAlleles), gene.name,gene.hugo,gene.strand))
                        self.referenceFile.write("EXONIC %d %s %s\n" % (len(exonic_locs),"/".join(transcript_groupings),"/".join([j for j in jxn_groups])))
                    else:
                        self.referenceFile.write("%s %s %s %s %s %s %s " % (self.chr,myLoc,refAllele,"/".join(altAlleles), gene.name,gene.hugo,gene.strand))
                        self.referenceFile.write("INTRONIC NA NA %s\n" % "-".join([str(s) for s in intronic_loc[3]]))


                    self.candFile.write(">%s_%s_%s_%s_%s_%s_%s_" % (self.chr,myLoc,gene.name,gene.hugo,gene.strand,"intronic","-".join([str(s) for s in intronic_loc[3]])))
                    self.candFile.write("REF-ALLELE_%s\n%s\n" %(refAllele,"".join([intronic_seq[0],refAllele,intronic_seq[-1]])))
                    n=0
                    for a in altAlleles:
                        n+=1
                        self.candFile.write(">%s_%s_%s_%s_%s_%s_%s_" % (self.chr,myLoc,gene.name,gene.hugo,gene.strand,"intronic","-".join([str(s) for s in intronic_loc[3]])))
                        self.candFile.write("ALT-ALLELE-%s_%s\n%s\n" %(str(n),a,"".join([intronic_seq[0],a,intronic_seq[-1]])))
                    k+=1
                    if k == len(snp_locs): return

            while k < len(snp_locs):

                myLoc = snp_locs[k]
                myKey = snp_key[myLoc]
                refAllele = myKey[0][0]
                altAlleles = myKey[0][1::]
                if myLoc + self.readLen >= len(self.seq): return
                intergenic_seq = ["".join(self.seq[(myLoc-1)-self.readLen:myLoc-1])]+[self.seq[myLoc-1]]+["".join(self.seq[myLoc:myLoc+self.readLen])]
                intergenic_loc = [["intergenic"], myLoc,"NULL",[myLoc-self.readLen,myLoc+self.readLen+1],intergenic_seq]
                self.referenceFile.write("%s %s %s %s NA NA NA " % (self.chr,myLoc,refAllele,"/".join(altAlleles)))
                self.referenceFile.write("INTERGENIC NA NA %s\n" % "-".join([str(s) for s in intergenic_loc[3]]))
                self.candFile.write(">%s_%s_NA_NA_NA_%s_%s_" % (self.chr,myLoc,"intergenic","-".join([str(s) for s in intergenic_loc[3]])))
                self.candFile.write("REF-ALLELE_%s\n%s\n" %(refAllele,"".join([intergenic_seq[0],refAllele,intergenic_seq[-1]])))
                n=0
                for a in altAlleles:
                    n+=1
                    self.candFile.write(">%s_%s_NA_NA_NA_%s_%s_" % (self.chr,myLoc,"intergenic","-".join([str(s) for s in intergenic_loc[3]])))
                    self.candFile.write("ALT-ALLELE-%s_%s\n%s\n" %(str(n),a,"".join([intergenic_seq[0],a,intergenic_seq[-1]])))
                k+=1











    def printCandidatesOnChromosome(self,TYPE='ALL'):
        geneTuples=[]

        for gene in self.genes:
            geneTuples.append((gene.start-1,gene.end))
            geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            geneInfo = gene.name+"|"+gene.hugo+"|"+gene.type+"|"+gene.chr+"|"+gene.strand+"|"
            if gene.end - gene.start > self.minLen:
                if self.findCands:
                    self.geneFile.write(">%s\n%s\n" % ( geneInfo+"GENE:"+str(gene.start)+"-"+str(gene.end), geneSeq))
                    for N in gene.novelJxns: self.candFile.write(">%s\n%s\n" % (geneInfo+"NJXN:"+listToString(N,["|","-"]),"".join([geneSeq[n[0]-gene.start:(n[1]-gene.start)+1] for n in N])))

            for start,end in gene.exons:    self.featureFile.write(">%s\n%s\n" % (geneInfo+"EXON:"+str(start)+"-"+str(end),geneSeq[start-gene.start:(end-gene.start)+1]))
            for I in gene.introns:   self.featureFile.write(">%s\n%s\n" % (geneInfo+"INTRON:"+listToString(I,["|","-"]),"".join([geneSeq[i[0]-gene.start:(i[1]-gene.start)+1] for i in I])))
            for J in gene.knownJxns: self.featureFile.write(">%s\n%s\n" % (geneInfo+"KJXN:"+listToString(J,["|","-"]),"".join([geneSeq[j[0]-gene.start:(j[1]-gene.start)+1] for j in J])))


        for g in geneTuples:
            if (g[1]-g[0]) > (self.readLen+1)*2:
                self.seq[g[0]+self.readLen:g[1]-self.readLen]=["N" for i in range((g[1]-self.readLen)-(g[0]+self.readLen))]

        chrLen=len(self.seq)


        self.headerFile.write("%s\n" % "\t".join(["@SQ","SN:"+self.chr,"LN:"+str(chrLen)]))
        k=0;BUFFER=200
        self.chrFile.write("%s\n" % (">"+self.chr))
        while True:
            self.chrFile.write("%s\n" % "".join([b for b in self.seq[BUFFER*k:BUFFER*(k+1)]]))
            k+=1
            if k*BUFFER > chrLen: break


############################################################################################################################################





############################################################################################################################################
####################################################### SIMULATION CODE ####################################################################
############################################################################################################################################

    def add_simulation_parameters(self,simulation_parameter_file):
        for line in open(simulation_parameter_file):
            line = line.split()
            if line[0] == 'cdna_length':         self.cdna_length  = int(line[-1])
            if line[0] == 'sonication_mean':     self.sonication_mean = int(line[-1])
            if line[0] == 'tail_adaptor'    :    self.tail_adaptor = line[-1]
            if line[0] == 'head_adaptor'    :    self.head_adaptor = line[-1]
            if line[0] == 'output_file_prefix':  self.prefix = line[-1]
            if line[0] == 'read_length'    :     self.readLen = int(line[-1])
            if line[0] == 'protocol'    :     self.protocol  = line[-1]
            if line[0] == 'size_selection_min'    :     self.size_selection_min  = int(line[-1])
            if line[0] == 'linear_amplification_values': self.linear_amplification_values = [int(x) for x in line[2::]]
            if line[0] == 'exponential_amplification_values': self.exponential_amplification_values = [int(x) for xi in line[2::]]




    def add_gene_key(self,gene_key_file):
        self.gene_key=dd(lambda: (1,0))
        for line in open(gene_key_file):
            line = line.split()
            self.gene_key[line[0]]=(int(line[-2]),int(line[-1]))











    def simulateGeneReadsOnChromosome(self,protocol="HUGO"):

        self.read_noise   = "".join(['A' for i in range(self.readLen)])
        self.qual = "".join(["I" for i in range(self.readLen)])
        self.readOutPut = sys.stdout

        for gene in self.genes:
            read_num = 0
            vals = self.gene_key[gene.hugo]
            gene.fragmentRNA(self.cdna_length,vals[0],vals[1])
            if self.protocol == "HUGO" or self.protocol == 'hugo':
                gene.linearAmplifyFrags(self.linear_amplification_values)
                gene.sonicateFrags(self.sonication_mean)

            if self.protocol == "OLEG" or self.protocol == 'oleg':
                gene.sonicateFrags(self.sonication_mean)
                gene.exponentialAmplifyFragments(self.exponential_amplification_values)



            geneSeq =  "".join([ base.capitalize() for base in self.seq[gene.start-1:gene.end]])
            for frag in gene.rnaFragments:
                if frag[-1] < self.size_selection_min: continue
                read_num+=1
                readID = ["@"+str(read_num),gene.chr,gene.name,frag[1],self.protocol]
                readSeq = "".join([geneSeq[x[0]-gene.start : (x[1]-gene.start)+1] for x in frag[2]])[0:self.readLen]
                read_coords= []
                read_dist,k   = 0,0
                while k < len(frag[2]):
                    if read_dist + (frag[2][k][1] - frag[2][k][0] ) + 1 <= len(readSeq):
                        read_coords.append(frag[2][k])
                        read_dist += (read_coords[-1][1] - read_coords[-1][0] ) + 1
                    else:
                        read_coords.append((frag[2][k][0],frag[2][k][0] + (len(readSeq) - read_dist)))
                        break
                    k+=1
                coord_strs = "G".join([str(x[0])+'-'+str(x[1]) for x in read_coords])
                readID.append(coord_strs)
                if len(readSeq) == self.readLen:    readID.append("N")
                else:
                    readID.append("Y")
                    if len(readSeq) < self.readLen: readSeq += self.tail_adaptor[0 : (self.readLen - len(readSeq))  ]
                    if len(readSeq) < self.readLen: readSeq = self.head_adaptor[max(0,len(self.head_adaptor)-(self.readLen-len(readSeq))) : len(self.head_adaptor) ]  + readSeq
                    if len(readSeq) < self.readLen: readSeq += self.read_noise[0: (self.readLen - len(readSeq)) ]

                for i in range(frag[0]):
                    #readID.append(str(i))
                    tmpID = "_".join(readID + [str(i)])
                    #readID_tmp = "_".join([j for j in readID])
                    self.readOutPut.write("%s\n%s\n%s\n%s\n" % (tmpID,readSeq,'+',self.qual))












