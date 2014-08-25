#!/usr/bin/env python



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
            self.geneID, self.tranID,   geneType, self.geneStatus, self.hugoName = line[9].split('"')[1], line[11].split('"')[1],line[13].split('"')[1], line[15].split('"')[1], line[17].split('"')[1]
            if geneType == "protein_coding" or geneType == "pc":
                self.geneType = "protein_coding"
            elif "pseudogene" in geneType:
                self.geneType = "pseudogene"
            elif "linc" in geneType:
                self.geneType = "lincRNA"
            elif "antisense" in geneType:
                self.geneType = "antisense"
            elif "IG_" == geneType[0:3] or "TR_" == geneType[0:3]:
                self.geneType = "immunoglobulin"
            elif geneType == "processed_transcript":
                self.geneType = "processed_transcript"
            elif "rRNA" in geneType:
                self.geneType = "ribosome"
            elif self.chr == "chrM" or self.chr == "MT" or self.chr == "M":
                self.geneType = "mitochondria"
            elif geneType == "miRNA":
                self.geneType = "miRNA"
            elif geneType == "snRNA" or geneType == "snoRNA":
                self.geneType = "snRNA"
            else:
                self.geneType = "misc"
            
         




