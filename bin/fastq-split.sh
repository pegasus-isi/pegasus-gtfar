#!/bin/bash

set -e
set -o pipefail

if [ $# -ne 2 ]; then
    echo "Usage: `basename $0` <READS> <SPLITS>"
    exit 1
fi

READS=$1
SPLITS=$2

if [ -z ${READS} ]; then
    echo "READS cannot be empty"
    exit 1
fi

if [ ! -e ${READS} ]; then
    echo "READS file does not exists"
    exit 2
fi

if [ -z ${SPLITS} ]; then
    echo "SPLITS cannot be empty"
    exit 1
fi


function split_reads {
    READS=$1
    SPLITS=$2

    myEXT=${READS##*.}
    FILENAME=$(echo ${READS} | awk -F\. '{print $1}')

    if [ ${myEXT} = gz ]; then

        myREADS=$(basename ${READS} ".${myEXT}")
        fEXT=${myREADS##*.}

        if [ ${fEXT} = fastq ] || [ ${fEXT} == fq ]; then

            gunzip --stdout ${READS} | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > ${READS}.tmp

        else

            echo "Need either a FastQ or FastQ.gz file"
            exit 3

        fi

    elif [ ${myEXT} == fastq ] || [ ${myEXT} == fq ]; then

        cat ${READS} | awk '{if (NR%4==0) printf $1"\n"; else printf $1" "}' > ${READS}.tmp

    else

        echo "READS file must have an extension of fastq, fq, fastq.gz, or fq.gz"
        exit 4

    fi

   LINE_CNT=$(wc -l ${READS}.tmp | awk '{print $1}')
   SPLIT_VAL=$(( $(( $LINE_CNT/$SPLITS )) + 1 ))
   split -d -l ${SPLIT_VAL} ${READS}.tmp
}

# Main
split_reads ${READS} ${SPLITS}
