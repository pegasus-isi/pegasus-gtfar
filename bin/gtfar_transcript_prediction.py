#!/usr/bin/env python


import os
import sys
from collections import defaultdict as dd


def process_file(INFO):

    exons = dd(lambda: dd(int))
    jxns  = dd(lambda: dd(int))
    for line in INFO:
        line = line.split()
        feature = line[0].split("@")[-1].split(":")[0]
        if feature not in ["EXON","KJXN"]: continue 
        gene = line[0].split("@")[0]
        coords = line[0].split("=")[-1]
        cnt = int(line[1])
        if feature == "EXON":
            exons[gene][coords]=cnt
        else:
            jxns[gene][coords]=cnt
    return exons,jxns
    

def reestimate_cov(latent_data,rate):
    estimates = dd(float)
    for group in latent_data:
        rate_sum = sum([rate[g] for g in group])
        cov = latent_data[group]
        print group,latent_data[group]
        for t in group:
            cov_estimates[t]+= cov * (rate[t] / (rate_sum) )
    return cov_estimates




def calculate_expected_reads_and_dist(my_group,my_rates,my_reads):
    my_shares = []
    for i in range(len(my_group)):
        my_shares.append(my_rates[i]/sum(my_rates))
    return [my_shares[i]*my_reads for i in range(len(my_shares))]
    


def EM_algorithm(latent_data,rate,dist):
    cov_estimates = dd(float)
    transcript_probs = dd(list)

    transcript_read_estimates = dd(float)
    print "INIT RATES"
    for r in rate:
        print r,rate[r]
    print ""
    for group in latent_data:
        transcript_shares=calculate_expected_reads_and_dist(group,[rate[g] for g in group],latent_data[group])
        for i in range(len(group)): transcript_read_estimates[group[i]]+=transcript_shares[i]
    print "NEW RATES"
    rate_dist = 0
    for r in rate:
        new_rate = transcript_read_estimates[r] / float(dist[r])
        rate_dist += (rate[r] - new_rate) * (rate[r] - new_rate)
        rate[r] = new_rate
        print r,rate[r]
    k=0
    while rate_dist >= 0.00001 and k<5:

        k+=1
        transcript_read_estimates = dd(float)
        for group in latent_data:
            transcript_shares=calculate_expected_reads_and_dist(group,[rate[g] for g in group],latent_data[group])
            for i in range(len(group)): transcript_read_estimates[group[i]]+=transcript_shares[i]
        print "ITER",k,"RATES"
        rate_dist = 0
        for r in rate:
            new_rate = transcript_read_estimates[r] / float(dist[r])
            rate_dist += (rate[r] - new_rate) * (rate[r] - new_rate)
            rate[r] = new_rate
            print r,rate[r]
    print rate_dist


    sys.exit()






def estimate_transcripts(gtf,rlen,exons,jxns):
    MODE = "INIT"
    for line in open(gtf):
        if line[0] == "#": continue
        line = line.split()
        if line[2] == "gene":
            if MODE == "EXON":
                exon_counts = exons[gene]
                jxn_counts = jxns[gene]
                if len(exon_counts) + len(jxn_counts) < 5 or len(jxn_counts)<3:
                    g=4
                    #print "NONE"
                else:

                    transcript_group_cnts,transcript_group_dists = dd(list),dd(list)
                    transcript_rates,transcript_dists = dd(float),dd(float)
                    for e in gene_exons:
                        group = tuple(sorted(gene_exons[e]))
                        dist = max(rlen,(float(int(e.split("-")[1])-int(e.split("-")[0]))) - rlen)
                        transcript_group_cnts[group].append(exon_counts[e])
                        transcript_group_dists[group].append(dist)
                        
                    multi_obs = []
                    for j in jxn_counts:
                        if len(j.split(","))>2:
                            transcript_union = gene_jxns[j.split("-")[0]]
                            for small_jxn in j.split("-")[1::]:
                                for t in transcript_union:
                                    if t not in gene_jxns[small_jxn]:
                                        transcript_union.remove(t)
                            group = tuple(sorted(transcript_union))
                            transcript_group_cnts[group].append(jxn_counts[j])
                            transcript_group_dists[group].append(rlen)
                            multi_obs.append(j.split("-")[0])
                            multi_obs.append(j.split("-")[1])
                            
                    for j in gene_jxns:
                        if j not in multi_obs:
                            group = tuple(sorted(gene_jxns[j]))
                            transcript_group_cnts[group].append(jxn_counts[j])
                            transcript_group_dists[group].append(rlen)
            
                    latent_data = {}
                    total_reads,total_span = 0,0
                    for g in transcript_group_cnts:
                        latent_data[g] = sum(transcript_group_cnts[g])
                        total_reads+= sum(transcript_group_cnts[g])
                        total_span += sum(transcript_group_dists[g])
                        for t in g:
                            transcript_dists[t]+=sum(transcript_group_dists[g])

                    print total_reads,total_span,total_reads / total_span 
                    ### FIX THIS SO THIS BECOMES THE PROBABILITY FOR EACH TRANSCRIPT ---> RELATIVE OF COURSE TO THIS NUMBER ###


                    for t in transcripts:          transcript_rates[t] = (total_reads/total_span) / len(transcripts)


                    EM_algorithm(latent_data,transcript_rates,transcript_dists)


            chrID,strand,gene = line[0],line[6],line[9].split('"')[1]+","+line[17].split('"')[1]
            MODE="GENE"
            transcripts = []
            gene_exons = dd(list)
            gene_jxns =   dd(list)
        elif line[2] == "transcript":
            tranID = line[11].split('"')[1]
            MODE = "TRANSCRIPT"
            transcripts.append(tranID)
        elif line[2] == "exon":
            if MODE == "TRANSCRIPT":
                coords = line[3]+'-'+line[4]
                gene_exons[coords].append(tranID)
                if strand == "-": last_loc = line[3]
                elif strand == "+": last_loc = line[4]
                MODE = "EXON"
               # if tranID == "ENST00000379319.1":
                #    print "YES start",coords,last_loc
            else:
                coords = line[3]+'-'+line[4]
                gene_exons[coords].append(tranID)
                if strand == "-":
                    jcoords = line[4]+","+last_loc
                    last_loc = line[3]
                elif strand == "+":
                    jcoords = last_loc+","+line[3]
                    last_loc = line[4]
                gene_jxns[jcoords].append(tranID)
               # if tranID == "ENST00000379319.1":
                #    print "YES exon",coords,last_loc
                 #   print "YES jcords",jcoords,last_loc

if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-g", "--gtf", default = None, type='string', help="GTF File Name")
    parser.add_option("-l", "--readlen", default = 100, type='int', help="GTF File Name")

    (options, args) = parser.parse_args()

    if not options.gtf:
        print "NEED GTF FILE"
        sys.exit()

    if len(args) == 1:
        FILE=open(args[0])
    else:
        FILE=sys.stdin

    exons,jxns = process_file(FILE)

    estimate_transcripts(options.gtf,options.readlen,exons,jxns)

