#!/bin/bash

### HARD CODED PATHS ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$PATH:$DIR



############################################################################################################################################################

NUMS=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40)

if [ $LENGTH -gt 101 ]; then
    TRIMS="50 75 100"
elif [ $LENGTH -gt 76 ]; then 
    TRIMS="50 75"
elif [ $LENGTH -gt 51 ]; then 
    TRIMS="50"
else
    TRIMS="0"
fi

TRIMSTR=$(echo $TRIMS | awk '{$1=$1}1' OFS=",")
TRIMARRAY=($TRIMS)

function count_and_split_reads {

    OUTDIR=$1; READS=$2; myPREFIX=$3
    CORES=$((  $(nproc) - 1 ))
    myEXT=${READS##*.}
    FILENAME=$(echo $READS | awk -F\. '{print $1}')

    if [ $OUTDIR == "NONE" ]; then OUTDIR="gtfar_$FILENAME"; fi 
        
    mkdir -p $OUTDIR; mkdir -p $OUTDIR"/filterDir" ; 
    if [ $myEXT = gz ]; then 
    	myREADS=$(basename $READS ".$myEXT"); fEXT=${myREADS##*.}
    	if [ $fEXT = fastq ] || [ $fEXT == fq ]; then
		    ln -fns $(readlink -m $READS) $OUTDIR/filterDir/READS.init 
		    cd $OUTDIR/filterDir   		
		    zcat READS.init | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > READS.tmp
        else 
		    error_quit "NEED FASTQ FILE OR FASTQ.GZ"; fi

    elif [ $myEXT == fastq ] || [ $myEXT == fq ]; then 
	    ln -fns $(readlink -m $READS) $OUTDIR/filterDir/READS.init 
	    cd $OUTDIR/filterDir   		
        cat READS.init | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > READS.tmp 
    else
        echo "Need read file extension to be fastq fq or fastq.gz or fq.gz"; exit; fi 
   
   LINE_CNT=$(wc -l READS.tmp | awk '{print $1}')
   SPLIT_VAL=$(( $(( $LINE_CNT/$CORES )) + 1 ))
   split -d -l $SPLIT_VAL READS.tmp

    ### READS ARE NOW SPLIT --- FILTER THEM NOW ###
    k=0
    printf "Filtering reads..."
    
    while [ $k -lt $CORES ]; do 
        pre_filter_fastq.py "x"${NUMS[$k]} -r $LENGTH -t $TRIMSTR -p reads"$k" & k=$(( $k + 1 ))
    done 
    wait
    if [ ! -s reads0.stats ]; then exit; fi
    cat reads*_reject.fastq > ../$myPREFIX.reject.fastq 
    cat reads*.stats | sort -k1n,1 | awk '{if ($1==L) {A+=$2;B+=$3} else {if (NR>1) print L,A,B;L=$1;A=$2;B=$3}} END{print L,A,B}' > ../$myPREFIX.adaptor.stats 
	printf "Success\n" 
    cd ../../
}






function error_quit {
    STATEMENT=$1
    echo $STATEMENT
    exit
}

function pass_fail {
    if [ $1 != 0 ]; then
        echo "FAIL"; exit
    fi
    if [ $# == 1 ]; then echo "SUCESS"; else echo -n $2; fi 
}




function token_cnt { echo $#
}

#############################################################################################################################################################
###############################################################  CHECKS  ####################################################################################
#############################################################################################################################################################













################################################################################################################################################################
################################################################################################################################################################

# ------------------------------------------ ALIGNMENT/PARSING  ------------------------------------------ #



     




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









