#!/bin/bash

### HARD CODED PATHS ###
export PATH=$PATH:/export/uec-gs1/knowles/analysis/tade/gtfar_source

#############################################################################################################################################################
##################################################### UTILITY FUNCTIONS #####################################################################################
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
    MYREADS=$READS; if [ -d $READS ]; then MYREADS=$READS/*.*; fi 
    for f in $MYREADS; do 
        FNAME=$(readlink -m $f);  FEXT="${FNAME##*.}"; FLINK=$WORKING/reads/$(basename $FNAME ".${FNAME##*.}")".fastq"
        if [ ! -f $FLINK ] && [ $(check_valid_extension $FNAME "READFILE") == "TRUE" ]; then
            ln -s $FNAME $FLINK 
            echo $PWD/$FLINK >> $WORKING/initial_reads.txt
        fi
    done 
    
    if [ ! -f $WORKING/initial_reads.txt ]; then
        echo "ERROR: NO VALID FASTQ FILES PROVIDES (extension must be fq or fastq)"; exit
    fi
    INIT_READS=initial_reads.txt 
    READLIST=$INIT_READS 
}


function load_parameters {
    PARAM_CLIPSTATUS="TRUE"; PARAM_CLIPLEN=$CLIPLENGTH; PARAM_READLEN=$LENGTH; PARAM_CLIPSEED=F1; PARAM_SUBS=$MISMATCHES
    if [ $LENGTH -lt 75 ]; then PARAM_CLIPSTATUS="FALSE"; fi 
    
    if [ $LENGTH -lt 65 ]; then 
        PARAM_INDEXLEN=$LENGTH
        PARAM_SEEDTYPE="F"$MISMATCHES
        if [ $MISMATCHES -gt 4 ]; then PARAM_SEEDTYPE="F4"; fi 
    elif [ $LENGTH -lt 128 ]; then 
        PARAM_INDEXLEN=$(( $LENGTH/2 )); PARAM_SEEDTYPE="F"$(( $MISMATCHES/2 ))
        if [ $(( $MISMATCHES%2 )) -gt 0 ]; then PARAM_SEEDTYPE="F"$(( $MISMATCHES/2+1 )); fi 
        if [ $MISMATCHES -gt 8 ]; then PARAM_SEEDTYPE="F8"; fi 
    else
        "echo READ LENGTH IS TOO LONG"; exit
    fi
    PARAM_CLIPSUBS=0
    if [ $MISMATCHES -gt 3 ]; then PARAM_CLIPSUBS=1; fi 
}







function get_seed {
    SUBS=$1;LEN=$2
    if   [ "$SUBS" -gt "6" ]; then echo "F4";
    elif [ "$SUBS" -gt "4" ]; then echo "F3"; 
    elif [ "$SUBS" -gt "2" ]; then echo "F2"; 
    elif [ "$SUBS" -gt "0" ]; then echo "F1";
    else echo "FO"
    fi
}













 
# 2) PREPARE READ DATA FROM A COMMAND LINE ARUGMENT  ------------------------------------------------------------------------------------------#

function check_valid_extension {
    TMP_FILE=$1; TMP_TYPE=$2; TMP_EXT="${TMP_FILE##*.}"
    if [ $TMP_TYPE == "READFILE" ]; then 
        if [ $TMP_EXT == "fastq" ] || [ $TMP_EXT == "fq" ]; then 
            echo "TRUE"
        else
            echo "FALSE"
        fi
    fi
}


# ---------------------------------------------------- PARAMETER SELECTION ---------------------------------------------------------------------- #




function get_ref_file {
    TMP_REF=$1; TMP_NAME=$(basename $1 ".${TMP_REF##*.}"); TMP_OUT=$TMP_REF
    if [ $2 != "GAPPED" ]; then 
        INDEX_SEARCH=$(dirname $TMP_REF)"/indexes/"$TMP_NAME"_"$PARAM_INDEXLEN"_"$PARAM_SEEDTYPE".index"
    else
        INDEX_SEARCH=$(dirname $TMP_REF)"/indexes/"$TMP_NAME"_"$PARAM_CLIPLEN"_"$PARAM_CLIPSEED".index"
    fi 
    if [ -f $INDEX_SEARCH ]; then echo $INDEX_SEARCH
    else                          echo $TMP_REF; fi
} 






    
function check_procs {
    READNUM=$(wc -l $1 | awk '{print $1}')
    if [ $READNUM -gt $(nproc) ]; then
        echo "WARNING: "$(nproc)" cores are being employed for" $READNUM "files; Speed may not be maximized"
    fi 
   }


function update_reads {
    REFNAME=$1
    NEWLIST=miss_"$REFNAME".txt
    if [ -f $NEWLIST ]; then rm $NEWLIST; fi 
    for i in $(cat $READLIST); do
        echo $PWD/reads/$(basename $i .fastq)"_miss_"$REFNAME".fastq" >> $NEWLIST            
    done
    READLIST=$NEWLIST
}


function check_mapping_inputs {
    TMPREF=$1; TMPREADS=$2
    if [ ! -f $TMPREF ]; then 
        echo "ERROR: "$TMPREF" does not exist"
        exit
    fi 
    for i in $(cat $READLIST); do 
        if [ ! -f $i ]; then 
            echo "ERROR: Invalid readfile "$i" input"
            exit
        fi
    done
}


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
















#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################


# ------------------------------------------ ALIGNMENT/PARSING  ------------------------------------------ #



function map_and_parse {
    REF=$1; TYPE=$2; TRUTH=$3; MYREF=$(get_ref_file $REF "UNGAPPED"); REFNAME=$(basename $REF ".${REF##*.}")
    
    check_mapping_inputs $MYREF $READLIST
    echo -n "RUNNING: perm $(basename $REF) $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME".log"....." 
    if [ "$TRUTH" != "NO" ]; then perm $MYREF $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME".log"  ; fi   
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    update_reads $REFNAME
    echo -n "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    for i in $(cat $INIT_READS); do 
        SAMPLE=$(basename $i .fastq)
        MYFILE=$REFNAME*_*$SAMPLE*.mapping
        MYPREF=$SAMPLE"_"$TYPE
        if [ "$TRUTH" != "NO" ]; then parse_alignment.py $MYFILE -p $MYPREF --strandRule $STRANDRULE & fi                                                    
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 200; done                                                     
    echo -n "."
    done
    wait 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    echo ""  
}
     

function novel_splice_search {
    REF=$1;TYPE=$2;TRUTH=$3; MYREF=$(get_ref_file $REF "GAPPED");  REFNAME=$(basename $REF ".${REF##*.}")
    
    check_mapping_inputs $MYREF $READLIST
    
    echo -n "RUNNING: clipR $(basename $REF) $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -e -s > $REFNAME"_clip.log.......""
    if [ "$TRUTH" != "NO" ]; then clipR $MYREF $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -e  --ignoreRepeatR 10 --ignoreDummyR 40 --noSamHeader -u -s > $REFNAME"_clip.log"; fi 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    
    echo -n "RUNNING: gtfar-splice-search {GAPFILES} > SPLICE_CANDIDATES.fa......"
    if [ "$TRUTH" != "NO" ]; then categorize_and_annotate_novel_splice.py $REFNAME*_*.sam -k $KEY --gtf $GTF -g $CHRS > SPLICE_CANDIDATES.fa ; fi 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    
    
    echo -n "RUNNING: perm SPLICE_CANDIDATES.fa $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > SPLICE_VALIDATE.log....." 
    if [ "$TRUTH" != "NO" ]; then perm SPLICE_CANDIDATES.fa $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u > SPLICE_VALIDATE.log  ; fi   
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    
    update_reads SPLICE_CANDIDATES
    echo -n "Running: gtfar-parse {MAPFILES} --strandRule $STRANDRULE...."
    for i in $(cat $INIT_READS); do 
        SAMPLE=$(basename $i .fastq)
        MYFILE=SPLICE_CANDIDATES*_*$SAMPLE*.mapping
        MYPREF=$SAMPLE"_"$TYPE
        if [ "$TRUTH" != "NO" ]; then parse_alignment.py $MYFILE -p $MYPREF  & fi                                                    
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 200; done                                                     
    echo -n "."
    done
    wait 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    echo ""  

}



function genomic_map {
    REF=$1; TYPE=$2; TRUTH=$3; MYREF=$(get_ref_file $REF "UNGAPPED"); MYCLIPREF=$(get_ref_file $REF "GAPPED"); REFNAME=$(basename $REF ".${REF##*.}")
    check_mapping_inputs $MYREF $READLIST
     
    echo -n "RUNNING: perm $(basename $REF) $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s > $REFNAME.log....." 
    if [ "$TRUTH" != "NO" ]; then perm $MYREF $READLIST --seed $PARAM_SEEDTYPE -v $PARAM_SUBS -B --printNM -u -s --outputFormat sam > $REFNAME.log  ; fi   
    if [ $? != 0 ]; then echo "FAILURE"; else echo "SUCESS"; fi  
    GENOMELIST=$READLIST
    update_reads $REFNAME
    

    check_mapping_inputs $MYREF $READLIST
    echo -n "RUNNING: clipR $(basename $REF) $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS -B --printNM -u -s > $REFNAME.log...." 
    if [ "$TRUTH" != "NO" ]; then clipR $MYCLIPREF $READLIST --seed $PARAM_CLIPSEED -v $PARAM_CLIPSUBS  --ignoreRepeatR 10 --ignoreDummyR 40 --noSamHeader -u -s > $REFNAME"_clip.log"; fi 
    if [ $? != 0 ]; then echo "FAILURE"; else echo "SUCESS"; fi  
    echo -n "Running: gtfar-parse {GENOMEFILES}...."
    for i in $(cat $INIT_READS); do 
        SAMPLE=$(basename $i .fastq)
        if [ "$TRUTH" != "NO" ]; then parse_genomic.py $REFNAME*$SAMPLE*.sam > $SAMPLE"_GENOME.stats" & fi 
        while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 200; done                                                     
        echo -n "."
    done 
    wait 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
    echo ""
}



#########################################################################################################################################################################
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


# ------------------------------------------ ANALYSIS/PARSING  ------------------------------------------ #


function combine_output {
    LASTREF=$1; TRUTH=$2; REFNAME=$(basename $REF ".${REF##*.}")
    echo -n "Running: gtfar-combineOutput {OUTPUTFILES}...."
    if [ "$TRUTH" != "NO" ]; then  
        for i in $(paste $GENOMELIST); do 
            SAMPLE=$(basename $i .fastq)
            PREV=$(echo $SAMPLE | awk -F\_miss_ '{print $1}')
            FIRSTMAP=$REFNAME*$SAMPLE.sam 
            SECONDMAP=$REFNAME*$SAMPLE*"_miss_"*.sam 
            if [ ! -f $FIRSTMAP ] || [ ! -f $SECONDMAP ]; then
                echo "ERROR: GENOME MAPPING WAS NOT COMPLETED"
                exit
            fi
            cat $FIRSTMAP $SECONDMAP > $PREV"_GENOME_vis.sam" & 
            cp reads/$SAMPLE"_miss_"$REFNAME"_miss_"$REFNAME* final_reads/$PREV"_unmapped.fastq" & 
            while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 50; done                                                     
            echo -n "."
        done
    fi
    wait

    if [ "$TRUTH" != "NO" ]; then  
        EXLOG=$(basename $EXONS ".${REF##*.}")".log"

        
        for i in $(cat $INIT_READS); do 
            SAMPLE=$(basename $i .fastq)
                combine_gtfar_files.py EXPRESSION  -e $SAMPLE"_EXONS_gene.cnts"  -i $SAMPLE"_INTRONS_gene.cnts"  -s $SAMPLE"_SPLICE_gene.cnts"  > results/$SAMPLE"_gene.cnts" & 
                combine_gtfar_files.py SPLICEDATA  -e $SAMPLE"_EXONS_splice.cnts"  -i $SAMPLE"_INTRONS_splice.cnts"  -s $SAMPLE"_SPLICE_splice.cnts"  > results/$SAMPLE"_splice.cnts" & 
                combine_gtfar_files.py STATS  -e $SAMPLE"_EXONS.stats"  -i $SAMPLE"_INTRONS.stats"  -s $SAMPLE"_SPLICE.stats" -g $SAMPLE"_GENOME.stats" -l $EXLOG  > results/$SAMPLE.stats & 
                #echo $SAMPLE"_GENOME_vis.sam" $SAMPLE"_EXONS_vis.sam" $SAMPLE"_INTRONS_vis.sam" $SAMPLE"_SPLICE_vis.sam"  results/$SAMPLE"_visualize.sam"
                cat $SAMPLE"_GENOME_vis.sam" $SAMPLE"_EXONS_vis.sam" $SAMPLE"_INTRONS_vis.sam" $SAMPLE"_SPLICE_vis.sam" > results/$SAMPLE"_visualize.sam" & 
                while [ $(jobs | wc -l) -gt $(nproc) ]; do sleep 200; done                                                     
                echo -n "."
        done 
    fi
    wait 
    if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "SUCESS"; fi     
}



function make_bam_files {
    ## FILES IN RESULTS DIR ##
    if [[ -z $(ls results/*visualize.sam 2> /dev/null) ]]; then 
        echo "Error: No Sam Files present to be converted to bam files"
    else
        echo "RUNNING: Production of Bam Files...."
        for sam in results/*_visualize.sam; do 
            PREF="${sam%.*}"
            samtools view -bS -o $PREF.bam $sam
            if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools sort $PREF.bam "$PREF"_sort
            if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            samtools index $PREF"_sort.bam"
            if [ $? != 0 ]; then echo "FAILURE"; exit; fi     
            mv $PREF.bam ./
        done
        if [ $? != 0 ]; then echo "FAILURE"; exit; else echo "Complete"; fi     
    fi
}









