#!/bin/bash

### HARD CODED PATHS ###
export PATH=$PATH:/export/uec-gs1/knowles/analysis/tade/gtfar_source



############################################################################################################################################################

function error_quit {
    STATEMENT=$1
    echo $STATEMENT
    exit
}

function pass_fail {
    if [ $1 != 0 ]; then
        echo "FAIL"; exit
    fi
    echo "SUCESS"
}




function count_tokens {
    echo $#
}

#############################################################################################################################################################
###############################################################  CHECKS  ####################################################################################
#############################################################################################################################################################


function check_valid_extension {
    TMP_FILE=$1; TMP_TYPE=$2; TMP_EXT="${TMP_FILE##*.}"
    if [ $TMP_TYPE == "READFILE" ]; then 
        if [ $TMP_EXT == "fastq" ] || [ $TMP_EXT == "fq" ]; then echo "TRUE"
        else                                                    echo "FALSE"; fi 
    fi
}

function check_input {
    if [ ! -f $1 ]; then
        error_quit "ERROR: File $1 not supplied"
    fi
}

function check_outputs {
    tmpREFNAME=$2; tmpREADLIST=$3; tmpTYPE=$4
    if [ $1 == "MAPPING" ]; then
        for t in $(cat $tmpREADLIST); do
            if [ ! -f $tmpREFNAME"_"$PARAM_PERMDATA"_"$(basename $t ".${t##*.}").mapping ] && [ ! -f $tmpREFNAME"_"$PARAM_PERMDATA"_"$(basename $t ".${t##*.}").sam ]; then echo "FALSE"; return; fi 
        done
    elif [ $1 == "MAPP" ]; then 
        echo "HI"
        echo $tmpREFNAME"_"$PARAM_PERMDATA"_"$tmpREADLIST
        echo "BYE"
        for t in $(cat $tmpREADLIST); do
            echo $tmpREFNAME"_"$PARAM_PERMDATA"_"$(basename $t ".${t##*.}")",mapping"
        done
        exit
    elif [ $1 == "SPLICEMAP" ]; then 
        for t in $(cat $tmpREADLIST); do
            PREF=$(basename $(basename $t ".${t##*.}") .genome)
            if [ ! -s $REFNAME"_"$PARAM_CLIPDATA"_"$PREF.sam ]; then echo "FALSE"; return; fi 
        done
    elif [ $1 == "FEATUREPARSE" ]; then
        for t in $(cat $tmpREADLIST); do
            tmpPREF=$(echo $(basename $t ".${t##*.}") | awk -F_miss '{print $1}')"_"$tmpTYPE
            if [ ! -s $tmpPREF".genecnts" ] || [ ! -s $tmpPREF".splicecnts" ] || [ ! -s $tmpPREF".vis" ] || [ ! -s $tmpPREF".stats" ]; then echo "FALSE"; return; fi 
        done
    elif [ $1 == "GENOMEPARSE" ]; then 
        for t in $(cat $tmpREADLIST); do 
            PREF=$(echo $(basename $i .fastq) | awk -F_miss '{print $1}')
            if [ ! -s $PREF"_GENOME.stats" ] || [ ! -s $PREF"_GENOME.vis" ]; then echo "FALSE"; return; fi
        done
    elif [ $1 == "COMBINE" ]; then 
        for t in $(cat $tmpREADLIST); do 
            PREF=$(basename $t .fastq)
            if [ ! -s results/$PREF".genecnts" ] || [ ! -s results/$PREF".splicecnts" ] || [ ! -s results/$PREF".stats" ] || [ ! -s results/$PREF"_visualize.sam" ]; then echo "FALSE"; return; fi 
        done
    fi    
    echo "TRUE"
}


function check_mapping_inputs {
    if [ $# != 2 ]; then error_quit "ERROR: Invalid Mapping Inputs $1"; fi 
    for i in $(cat $2); do 
        if [ ! -f $i ]; then error_quit "ERROR: Invalid Mapping Inputs"; fi  
    done
}


function check_procs {
    READNUM=$(wc -l $1 | awk '{print $1}')
    if [ $READNUM -gt $(nproc) ]; then  echo "WARNING: "$(nproc)" cores are being employed for" $READNUM "files; Speed may not be maximized"; fi 
    if [ $(( READNUM - 5 )) -gt $(nproc) ]; then  error_quit "ERROR: Far too many read files for number of supplied cores $READNUM vs $(nproc)"; fi
   }

#############################################################################################################################################################
##################################################### DATA LOADING      #####################################################################################
#############################################################################################################################################################


# 1) LOAD ANNOTATION/SPECIES DATA FROM A COMMAND LINE ARUGMENT  ------------------------------------------------------------------------------------------#
function load_annotation {
    if [ $ANNOTATION == "gc18" ] || [ $ANNOTATION == "GC18" ]; then
        SPECIES="HUMAN"; DATA=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human; TRANSCRIPTOME=$DATA/gencode18/readLength"$LENGTH"
        EXONS=$TRANSCRIPTOME/*exonSeqs.fa; INTRONS=$TRANSCRIPTOME/*intronSeqs.fa ; GENES=$TRANSCRIPTOME/*geneSeqs.fa ; KEY=$TRANSCRIPTOME/*.key 
        GTF=$DATA/gencode18/gencode.v18.annotation.gtf; CHRS=$DATA/GRCh37_ensemble19/chrs; GENOME=$DATA/GRCh37_ensemble19/GRCh37.p11.genome.fa 

    elif [ $ANNOTATION == "MONKEY" ] || [ $ANNOTATION == "monkey" ]; then
        SPECIES="MONKEY"; DATA=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/monkey; TRANSCRIPTOME=$DATA/norgren_transcriptome/readLength"$LENGTH"
        EXONS=$TRANSCRIPTOME/*exonSeqs.fa; INTRONS=$TRANSCRIPTOME/*intronSeqs.fa ; GENES=$TRANSCRIPTOME/*geneSeqs.fa ; KEY=$TRANSCRIPTOME/*.key 
        GTF=$DATA/norgren_transcriptome/norgren_convert.gtf; CHRS=$DATA/norgren_genome/chrs; GENOME=$DATA/norgren_genome/norgren_genome.txt  
    else
        show_help "ANNOTATION $ANNOTATION NOT SUPPORTED -choose gc17 or nothing"
    fi
}
# 1) END --------------------------------------------------------------------------------------------------------------------------------------------------#

function load_reads {
    if [ ! -d $WORKING/reads ]; then  mkdir $WORKING/reads; fi     
    MYREADS=$READS; if [ $(count_tokens $READS) == 1 ] && [ -d $READS ]; then MYREADS=$READS/*.*; fi  
    
    for f in $MYREADS; do 
        FNAME=$(readlink -m $f);  FEXT="${FNAME##*.}"; FLINK=$WORKING/reads/$(basename $FNAME ".${FNAME##*.}")".fastq"
        if [ ! -f $FLINK ] && [ $(check_valid_extension $FNAME "READFILE") == "TRUE" ]; then
            ln -s $FNAME $FLINK; echo $PWD/$FLINK >> $WORKING/initial_reads.txt
        fi
    done     
    if [ ! -f $WORKING/initial_reads.txt ]; then error_quit "ERROR: No valid fastq read files provided (extension must be fq or fastq)"; fi 
    INIT_READS=initial_reads.txt 
    READLIST=$INIT_READS 
}

function load_parameters {
    PARAM_CLIPSTATUS="TRUE"; PARAM_CLIPLEN=$CLIPLENGTH; PARAM_READLEN=$LENGTH; PARAM_CLIPSEED=F1; PARAM_SUBS=$MISMATCHES
    if [ $LENGTH -lt 75 ]; then PARAM_CLIPSTATUS="FALSE"; fi 
    if [ $LENGTH -lt 65 ]; then 
        PARAM_INDEXLEN=$LENGTH; PARAM_SEEDNUM=$MISMATCHES; PARAM_SEEDTYPE="F"$MISMATCHES; if [ $MISMATCHES -gt 4 ]; then PARAM_SEEDTYPE="F4"; fi 
    elif [ $LENGTH -lt 128 ]; then 
        PARAM_INDEXLEN=$(( $LENGTH/2 ));    PARAM_SEEDTYPE="F"$(( $MISMATCHES/2 ))
        if [ $(( $MISMATCHES%2 )) -gt 0 ];  then  PARAM_SEEDNUM=$(( MISMATCHES/2+1 )); PARAM_SEEDTYPE="F"$PARAM_SEEDNUM; fi 
        if [ $MISMATCHES -gt 8          ];  then  PARAM_SEEDTYPE="F4"; fi 
    else
        error_quit "ERROR: Read length is too long"
    fi
    PARAM_CLIPSUBS=0; if [ $MISMATCHES -gt 3 ]; then PARAM_CLIPSUBS=1; fi 
    PARAM_PERMDATA="B_"$PARAM_SEEDNUM"_"$MISMATCHES; PARAM_CLIPDATA="A_1_"$PARAM_CLIPSUBS"_"$PARAM_CLIPLEN
}


function load_ref_file {
    TMP_REF=$1; TMP_NAME=$(basename $1 ".${TMP_REF##*.}"); TMP_OUT=$TMP_REF
    if [ $2 != "GAPPED" ]; then 
        INDEX_SEARCH=$(dirname $TMP_REF)"/indexes/"$TMP_NAME"_"$PARAM_INDEXLEN"_"$PARAM_SEEDTYPE".index"
    else
        INDEX_SEARCH=$(dirname $TMP_REF)"/indexes/"$TMP_NAME"_"$PARAM_CLIPLEN"_"$PARAM_CLIPSEED".index"
    fi 
    if [ -f $INDEX_SEARCH ]; then echo $INDEX_SEARCH
    else                          echo $TMP_REF; fi
} 







 
# 2) PREPARE READ DATA FROM A COMMAND LINE ARUGMENT  ------------------------------------------------------------------------------------------#



# ---------------------------------------------------- PARAMETER SELECTION ---------------------------------------------------------------------- #





function update_reads {
    REFNAME=$1 
    for i in $(cat $READLIST); do
        echo $PWD/reads/$(basename $i .fastq)"_miss_"$REFNAME".fastq" >> NEW_READS.tmp         
    done
    mv NEW_READS.tmp miss_"$REFNAME".txt 
    READLIST=miss_"$REFNAME".txt
}




################################################################################################################################################################
################################################################################################################################################################

# ------------------------------------------ ALIGNMENT/PARSING  ------------------------------------------ #



     
function perm_map {
    REF=$1; TYPE=$2; MYREF=$(load_ref_file $REF "UNGAPPED"); REFNAME=$(basename $REF ".${REF##*.}")
    echo -n "RUNNING: perm $(basename $REF) $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME".log"....." 
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "MAPPING" $REFNAME $READLIST) == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else perm $MYREF $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u  > $REFNAME".log"; pass_fail $? 
    fi
}




function feature_parse {
    REF=$1; TYPE=$2; TRUTH=$3; REFNAME=$(basename $REF ".${REF##*.}")
    echo -n "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "FEATUREPARSE" $REFNAME $READLIST "EXONS") == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else
        for i in $(cat $READLIST); do 
            READNAME=$(basename $i ".${i##*.}"); MAPDATA="B_"$PARAM_SEEDNUM"_"$MISMATCHES; MAPFILE=$REFNAME"_"$MAPDATA"_"$READNAME".mapping"; PREF=$(echo $READNAME | awk -F_miss '{print $1}')
            check_input $MAPFILE; echo -n "."
            parse_alignment.py $MAPFILE -p $PREF"_"$TYPE --strandRule $STRANDRULE &
            while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
        done
        wait; pass_fail $?
    fi
    update_reads $REFNAME
}


function splice_find_and_conditional_map {
    REF=$1;TYPE=$2; MYREF=$(load_ref_file $REF "GAPPED");  REFNAME=$(basename $REF ".${REF##*.}")

    if [ $LENGTH -gt 74 ]; then 
        check_mapping_inputs $MYREF $READLIST    
        echo -n "RUNNING: clipR $(basename $REF) $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -e -s > $REFNAME"_clip.log.......""
        if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "SPLICEMAP" $REFNAME $READLIST) == "TRUE" ]; then echo "SKIPPING (previously completed)"
        else clipR $MYREF $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -e  --ignoreRepeatR 10 --ignoreDummyR 40 --noSamHeader > $REFNAME"_clip.log"; pass_fail $?; fi 
        echo -n "RUNNING: gtfar-splice-search {GAPFILES} > SPLICE_CANDIDATES.fa......"
        if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ -f SPLICE_CANDIDATES.fa ]; then echo "SKIPPING (previously completed)"
        else categorize_and_annotate_novel_splice.py $REFNAME*_*.sam -k $KEY --gtf $GTF -g $CHRS > SPLICE_CANDIDATES.fa ;  pass_fail $?; fi 
    fi 
    
    #echo $READLIST
    if [ ! -s SPLICE_CANDIDATES.fa ]; then
        echo -n "RUNNING: perm SPLICE_CANIDATES.fa $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME.log....." 
        echo "SKIPPING (no splice candidates found)"
        for i in $(cat $INIT_READS); do PREF=$(basename $i .fastq); touch $PREF"_SPLICE.stats"; touch $PREF"_SPLICE.genecnts"; touch $PREF"_SPLICE.splicecnts"; touch $PREF"_SPLICE.vis"; done
    elif [ ! -z $3 ] && [ $3 == "SKIP" ] &&  [ $(check_outputs "MAPPING" SPLICE_CANDIDATES $READLIST) == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else    perm_map SPLICE_CANDIDATES.fa "SPLICE"; feature_parse SPLICE_CANDIDATES.fa "SPLICE"; fi 
} 





function perm_genome_map {
    REF=$1; TYPE=$2; MYREF=$(load_ref_file $REF "UNGAPPED"); MYCLIPREF=$(load_ref_file $REF "GAPPED"); REFNAME=$(basename $REF ".${REF##*.}")
    check_mapping_inputs $MYREF $READLIST 
    echo -n "RUNNING: perm $(basename $REF) $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME.log....." 
    
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "MAPPING" $REFNAME $READLIST) == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else perm $MYREF $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s --outputFormat sam > $REFNAME.log; pass_fail $?; fi
    GENOMELIST=$READLIST; update_reads $REFNAME; 
    check_mapping_inputs $MYCLIPREF $READLIST

    echo -n "RUNNING: clipR $(basename $REF) $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -B --printNM -u -s > $REFNAME.log...." 
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "SPLICEMAP" $REFNAME $READLIST) == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else clipR $MYCLIPREF $READLIST --anchorL $PARAM_CLIPLEN --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS  --ignoreRepeatR 10 --ignoreDummyR 40 --noSamHeader -u -s > $REFNAME".cliplog"; pass_fail $?; fi
    update_reads $REFNAME 
}


function genome_parse {
    REF=$1; TYPE=$2;  MYREF=$(load_ref_file $REF "UNGAPPED"); MYCLIPREF=$(load_ref_file $REF "GAPPED"); REFNAME=$(basename $REF ".${REF##*.}")
    echo -n "Running: gtfar-parse {GENOMEFILES}...."
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "GENOMEPARSE" $REFNAME $GENOMELIST "GENOME") == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else
        for i in $(cat $GENOMELIST); do
            MIDNAME=$(basename $i ".${i##*.}");  PREF=$(echo $(basename $i) | awk -F_miss '{print $1}')
            F1=$REFNAME"_"$PARAM_PERMDATA"_"$MIDNAME.sam; F2=$REFNAME"_"$PARAM_CLIPDATA"_"$MIDNAME"_miss_"*sam 
            if [ $(count_tokens $F1 $F2) != 2 ]; then error_quit "WRONG FILE NUMBER"; fi 
            parse_genomic.py -p $PREF"_"$TYPE $F1 $F2 & 
        done
        wait
        pass_fail $?
    fi 
}



#########################################################################################################################################################################
# ------------------------------------------ COMBINING  ------------------------------------------ #



function combine_output {
    REF=$1; TMPREADS=$2; REFNAME=$(basename $REF ".${REF##*.}"); EXLOG=$(basename $EXONS ".${REF##*.}")".log"
    echo -n "Running: gtfar-combineOutput {OUTPUTFILES}...."
    if [ ! -z $3 ] && [ $3 == "SKIP" ] && [ $(check_outputs "COMBINE" $REFNAME $INIT_READS) == "TRUE" ]; then echo "SKIPPING (previously completed)"
    else
        for i in $(cat $TMPREADS); do 
            SAMPLE=$(basename $i .fastq) 
            cp $i final_reads/$SAMPLE"_unmapped".fastq & 
            combine_gtfar_files.py EXPRESSION  -e $SAMPLE"_EXONS.genecnts"  -i $SAMPLE"_INTRONS.genecnts"  -s $SAMPLE"_SPLICE.genecnts"  > results/$SAMPLE".genecnts" & 
            combine_gtfar_files.py SPLICEDATA  -e $SAMPLE"_EXONS.splicecnts"  -i $SAMPLE"_INTRONS.splicecnts"  -s $SAMPLE"_SPLICE.splicecnts"  > results/$SAMPLE".splicecnts" & 
            combine_gtfar_files.py STATS  -e $SAMPLE"_EXONS.stats"  -i $SAMPLE"_INTRONS.stats"  -s $SAMPLE"_SPLICE.stats" -g $SAMPLE"_GENOME.stats" -l $EXLOG  > results/$SAMPLE.stats & 
            cat $SAMPLE"_GENOME.vis" $SAMPLE"_EXONS.vis" $SAMPLE"_INTRONS.vis" $SAMPLE"_SPLICE.vis" > results/$SAMPLE"_visualize.sam" & 
            while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 10; done                                                     
            echo -n "."
        done
        wait
        pass_fail $? 
    fi
}

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


function check_dir {
    DIR=$1
    if [ -d $DIR ]; then echo "WARNING: working directory ("$DIR") already exists; Files may be overwritten"; 
    else mkdir $DIR; fi 
}

function make_directories {
    if [ -d $WORKING ]; then echo "WARNING: working directory ("$WORKING") already exists; Files may be overwritten"; 
    else mkdir $WORKING; fi 
    if [ -d $OUTPUT ]; then echo "WARNING: output directory ("$OUTPUT") already exists; Files may be overwritten"; 
    else mkdir $OUTPUT; fi 
    if [ ! -d $WORKING/results ]; then  mkdir $WORKING/results; fi 
    if [ ! -d $WORKING/final_reads ]; then  mkdir $WORKING/final_reads; fi 
}

function start_up_output {
    echo ""   
    echo "GTFAR RNA Pipeline Verson 0.0.5"
    echo "SPECIES (-s):            "$SPECIES
    echo "GTF ANNOTATION (-g):     "$ANNOTATION "( "$(basename $GTF) ")"
    echo "Reads:                   "$READS
    echo "Substitutions allowed:   "$MISMATCHES         
    echo "Working Directory:       "$WORKING
    echo "Output Directory:        "$OUTPUT
    echo "Number of readFiles:     "$(cat $READLIST | wc -l)
    echo "Number of processors:    "$(nproc)
    echo ""
}











