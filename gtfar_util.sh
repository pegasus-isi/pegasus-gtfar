#!/bin/bash

### HARD CODED PATHS ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$PATH:$DIR



############################################################################################################################################################








function gtfar_fail {
    STATEMENT="$1"
    printf "$STATEMENT\n"
    exit
}

function gtfar_check_for_errors {
    if [ $# -gt 1 ]; then STATEMENT="$2"; else STATEMENT="FAILURE"; fi 
    if [ $1 != 0 ];  then printf $STATEMENT"\n"; fi 
}



function check_progress {
    if [ $1 != 0 ]; then echo "FAILURE"; exit; fi
}




function token_cnt { echo $#
}

#############################################################################################################################################################
###############################################################  CHECKS  ####################################################################################
#############################################################################################################################################################













################################################################################################################################################################
################################################################################################################################################################

# ------------------------------------------ ALIGNMENT/PARSING  ------------------------------------------ #



     
function perm_map {
    
    myREF=$1; myTYPE=$2; EXT="${myREF##*.}"
    echo -n "RUNNING: perm $(basename $myREF) $(basename $READLIST) --seed $SEED -v $MISMATCHES -B --printNM -u -s > $myTYPE.log....." # REFNAME".log"....." 
  

    if [ -f $myTYPE.log ] && [ $SKIP == "TRUE" ]; then echo "SKIPPING (previously completed)";
    else
        if [ $myTYPE != "AGENOME" ]; then  
            if [ $EXT == "index" ]; then  printf "...index found..."; perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM -u > $myTYPE.log; pass_fail $?
            else  perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"$SEED".index > $myTYPE.log; pass_fail $?; fi 
        else
            if [ $EXT == "index" ]; then  printf "...index found..."; perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM --outputFormat sam -u > $myTYPE.log; pass_fail $?
            else  perm $myREF $READLIST --seed $SEED -v $MISMATCHES -B --printNM --outputFormat sam -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"$SEED".index > $myTYPE.log; pass_fail $?; fi 
        fi
    fi
}

function clip_map {
    myREF=$1; myTYPE=$2; EXT="${myREF##*.}"
    echo -n "RUNNING: clipR $(basename $myREF) $(basename $READLIST) --seed F1 -v $MISMATCHES -B --printNM -u -s > $myTYPE.log....." # REFNAME".log"....." 
   
    
    if [ $LENGTH -lt 75 ]; then echo "SKIPPING (readlength is less than 75bp)"
    elif [ -f $myTYPE.log ] && [ $SKIP == "TRUE" ]; then printf "SKIPPING (previously completed)\n"; 
    elif [ $EXT == "index" ]; then printf "...index found...";  clipR $myREF $READLIST --seed F1 -v 1 -e -u -s --ignoreDummyR 40 --ignoreRepeatR 10 --noSamHeader > $myTYPE.log; pass_fail $?; 
    else clipR $myREF $READLIST --seed F1 -v 1 -e -u -s refs/gtfar-"$LENGTH"-"$myTYPE"-"F1.index" --ignoreDummyR 40 --ignoreRepeatR 10 --noSamHeader > $myTYPE.log; pass_fail $?; fi    
}

function perm_parse {
    myREF=$1;  RNAME=$(basename $myREF ".${myREF##*.}"); MAP_PREFIX=$RNAME"_B_"$SEEDSCR"_"$MISMATCHES; myTYPE=$2
    
    printf  "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    rm -f reads/TMP_READS.txt 
    for i in $(cat $READLIST); do 
        #MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(basename $(basename $i "_miss_${i##*_miss_}") .fastq); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq
        MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(echo $(basename $i .fastq) | awk -F_ '{print $1}'); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq

        if [ -f $PARSE_PREFIX"_"$myTYPE.vis ] && [ $SKIP == "TRUE" ]; then echo -n "SKIPPING,"; 
        elif [ $(token_cnt $MAPFILE) -eq 1 ]; then parse_alignment.py $MAPFILE --strandRule $STRANDRULE > $PARSE_PREFIX"_"$myTYPE.vis & 
        else error_quit "ERROR: AMIBIGUOUS COPIES OF MAPPING FILES"; fi  
        
        if [ -f $NEW_READS ]; then echo $NEW_READS >> reads/TMP_READS.txt; fi  
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
    done
    wait; pass_fail $?
    if [ -f reads/TMP_READS.txt ]; then mv reads/TMP_READS.txt reads/readlist.txt; fi 
}

function clip_parse {
    myREF=$1; myTYPE=$2; RNAME=$(basename $myREF ".${myREF##*.}"); MAP_PREFIX=$RNAME"_A_1_1_"$CLIPLEN
    printf  "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    rm -f TMP_READS.txt 
    for i in $(cat $READLIST); do 
        MAPFILE=$MAP_PREFIX"_"$(basename $i ".${i##*.}").*[mapping,sam]; PARSE_PREFIX=$(echo $(basename $i) | awk -F_ '{print $1}'); NEW_READS=$OUT/reads/$(basename $i .fastq)"_miss_"$RNAME.fastq
        
        if [ -f $PARSE_PREFIX"_"$myTYPE.vis ] && [ $SKIP == "TRUE" ]; then echo -n "SKIPPING,"; 
        elif [ $(token_cnt $MAPFILE) -eq 1 ]; then parse_alignment.py $MAPFILE --strandRule $STRANDRULE > $PARSE_PREFIX"_"$myTYPE.vis & 
        else error_quit "ERROR: AMIBIGUOUS COPIES OF MAPPING FILES"; fi  
        
        if [ -f $NEW_READS ]; then echo $NEW_READS >> reads/TMP_READS.txt; fi  
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
    done
    wait; pass_fail $?
    if [ -f reads/TMP_READS.txt ]; then mv reads/TMP_READS.txt reads/readlist.txt; fi 
}



function combine_vis_files {

for i in $(cat reads/read_index.table); do 
    RNAME=$(echo $i | awk -F\=== '{print $1}')
    RFULL=$(echo $i | awk -F\=== '{print $2}')
    cat $RNAME*.vis > $RFULL.sam 
done
}


#########################################################################################################################################################################
# ------------------------------------------ COMBINING  ------------------------------------------ #




function make_bam_files {
    if [[ -z $(ls results/*visualize.sam 2> /dev/null) ]]; then 
        echo "Error: No Sam Files present to be converted to bam files"
    else
        echo "RUNNING: Production of Bam Files...."
        for sam in results/*_visualize.sam; do 
            PREF="${sam%.*}"
            samtools view -bS -o $PREF.bam $sam;  if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools sort $PREF.bam "$PREF"_sort; if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools index $PREF"_sort.bam";      if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            mv $PREF.bam ./
        done
        if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "Complete"; fi     
    fi
}


#########################################################################################################################################################################









