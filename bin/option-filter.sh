#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: option-filter.sh <PREFIX> <READS> <LENGTH>"
    exit 1
fi

PREFIX=$1
READS=$2
LENGTH=$3

if [ -z $PREFIX ]; then
    echo "PREFIX cannot be empty"
    exit 1
fi

if [ -z $READS ]; then
    echo "READS cannot be empty"
    exit 1
fi

if [ ! -e $READS ]; then
    echo "READS file does not exists"
    exit 2
fi

if [ -z $LENGTH ]; then
    echo "LENGTH cannot be empty"
    exit 1
fi


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

function count_and_split_reads {
    myPREFIX=$1
    READS=$2

    CORES=$((  $(nproc) - 1 ))
    myEXT=${READS##*.}
    FILENAME=$(echo $READS | awk -F\. '{print $1}')

    mkdir filterDir

    if [ $myEXT = gz ]; then
        myREADS=$(basename $READS ".$myEXT"); fEXT=${myREADS##*.}
        if [ $fEXT = fastq ] || [ $fEXT == fq ]; then
            ln -fns $(readlink -m $READS) filterDir/READS.init
            cd filterDir
            zcat READS.init | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > READS.tmp
        else
            echo "NEED FASTQ FILE OR FASTQ.GZ"
             exit 3
        fi

    elif [ $myEXT == fastq ] || [ $myEXT == fq ]; then
        ln -fns $(readlink -m $READS) filterDir/READS.init
        cd filterDir
        cat READS.init | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > READS.tmp

    else
        echo "READS file must have an extension of fastq, fq, fastq.gz, or fq.gz"
        exit 4
    fi

   LINE_CNT=$(wc -l READS.tmp | awk '{print $1}')
   SPLIT_VAL=$(( $(( $LINE_CNT/$CORES )) + 1 ))
   split -d -l $SPLIT_VAL READS.tmp

    ### READS ARE NOW SPLIT --- FILTER THEM NOW ###
    k=0
    echo "Filtering reads..."

    while [ $k -lt $CORES ]; do
        pre_filter_fastq.py "x"${NUMS[$k]} -r $LENGTH -t $TRIMSTR -p reads"$k" &
        k=$(( $k + 1 ))
    done
    wait

    if [ ! -s reads0.stats ]; then
        exit
    fi
    cat reads*_reject.fastq > ../$myPREFIX.reject.fastq
    cat reads*.stats | sort -k1n,1 | awk '{if ($1==L) {A+=$2;B+=$3} else {if (NR>1) print L,A,B;L=$1;A=$2;B=$3}} END{print L,A,B}' > ../$myPREFIX.adaptor.stats
    echo "Success"
    cd ..
}

# Main
count_and_split_reads $PREFIX $READS
