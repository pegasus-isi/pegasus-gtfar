#!/bin/bash

set -e
set -o pipefail

if [ $# -lt 2 ]; then
    echo "Usage: `basename $0` <REGEX> <OUTPUT>"
    exit 1
fi

REGEX=${@:1:$#}
OUTPUT=${!#}

if [ -z "${REGEX}" ]; then
    echo "REGEX cannot be empty"
    exit 2
fi

if [ -z "${OUTPUT}" ]; then
    echo "OUTPUT cannot be empty"
    exit 3
fi

cat ${REGEX} | sort -k1n,1 | awk '{if ($1==L) {A+=$2;B+=$3} else {if (NR>1) print L,A,B;L=$1;A=$2;B=$3}} END{print L,A,B}' > ${OUTPUT}

exit $?
