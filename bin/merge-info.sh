#!/bin/bash

set -e
set -o pipefail

if [ $# -lt 2 ]; then
    echo "Usage: `basename $0` <REGEX> <OUTPUT>"
    exit 1
fi

REGEX=${@:0:$#}
OUTPUT=${!#}

if [ -z "${REGEX}" ]; then
    echo "REGEX cannot be empty"
    exit 2
fi

if [ -z "${OUTPUT}" ]; then
    echo "OUTPUT cannot be empty"
    exit 3
fi

cat ${REGEX} | sort | awk '{if ($1==C && $3==T && $4==S && $5==E) {n=substr($16,2,3);D[16]+=n;D[22]+=n} else {printf D[1]; for (i=2;i<=NF;i++) printf "\t"D[i]; printf "\n"; for (i=1;i<=NF;i++) D[i]=$i; n=substr($16,2,3)+0; D[16]=n;D[22]=n;C=$1;T=$3;S=$4;E=$5}} END{printf D[1]; for (i=2;i<=25;i++) printf "\t"D[i]; printf "\n"}' > ${OUTPUT}

exit $?
