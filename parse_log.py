#!/usr/bin/env python


import sys
import os
				
def process_file(logfile,reads,OUTSTR):

    INIT=True; k=0; TOT=0
    f=open(logfile)
    subs=[]
    for line in f:
        line=line.split()
        
        if len(line) < 3: continue 

        if line[1] == 'effective': RL = line[5]
        elif line[3] == 'allowed':   SUBS = line[0]

        elif line[2] == "Reference":
            REF=line[3].split("/")[-1]

        elif line[0] == "Reference":
            REF=line[1].split("/")[-1]

        elif line[1] == 'Reads:,':
            READS=int(line[5].strip(',')); MAPPED=int(line[7].strip(','))
    

            if INIT:
                INIT=False
                #print READS,MAPPED
                START_READS=READS
                TOT+=MAPPED
                print "INFO: STARTING READS: ",READS," ( length=",RL,")"
                print "INFO: Round",1," | MAPPED",MAPPED,"%",float(MAPPED)/START_READS
            else:
                k+=1
                TOT+=MAPPED
                print "INFO: Round",k," | MAPPED",MAPPED,"%",float(TOT)/START_READS


    #print READS,MAPPED,subs



if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)

    parser.add_option("-z", "--hg19", action = 'store_true', default = False, help="Gapped input")
    parser.add_option("-s", "--strand", default = None, type='string', help="0,+, OR 16,-")
    parser.add_option("-o", "--outstring", default = "YO", type='string', help="0,+, OR 16,-")
    parser.add_option("-r", "--reads", default = 0, type=int, help="inital reads")

    (options, args) = parser.parse_args()

    if len(args)==0:
        print "PARSER FOR LOG FILE (FRESH FILE)"
        print ""
    elif len(args)==1:
        logFile=args[0]
        process_file(logFile,options.reads,options.outstring)
    else:
        print "TOO MANY ARGS"
        sys.exit()
