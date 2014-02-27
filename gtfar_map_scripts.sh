#!/bin/bash

### HARD CODED PATHS ###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$PATH:$DIR
source $DIR/gtfar_util.sh



############################################################################################################################################################


function map_and_parse_reads { 
    MAPDIR=$1; READS=$2; TAG=$3
    mkdir -p $MAPDIR 
   

    setup_perm_seeds
    setup_perm_reads
    setup_perm_refs
    SKIP="TRUE"
    SKIP="FALSE"
    cd $MAPDIR 
    map_and_parse_reads_to_features 
    map_and_parse_reads_to_genome
    if [ -z $SPLICES ]; then
        read_cands=(reads*miss_FEATURES_miss_GENOME.fastq)
        if [ -e $read_cands ]; then for r in reads*miss_FEATURES_miss_GENOME.fastq; do mv $r $(basename $r | awk -F\_miss '{print $1}')"_unmapped.fastq"; done; fi
    else    
        map_and_parse_reads_to_splices 
        read_cands=(reads*miss_FEATURES_miss_GENOME_miss_SPLICES.fastq)
        if [ -e $read_cands ]; then for r in reads*miss_FEATURES_miss_GENOME_miss_SPLICES.fastq; do mv $r $(basename $r | awk -F\_miss '{print $1}')"_unmapped.fastq"; done ; fi
    fi
    cd ..
}

function clip_and_parse_reads {
    CLIPDIR=$1; READS=$2; SEED="F1"
    setup_clip_seed
    setup_clip_reads $CLIPDIR $READS
    setup_clip_refs
    cd $CLIPDIR
    clipR $myGENES gene_reads.txt --seed $SEED --anchorL $anchor -e -v $MISMATCHES -s -u --noSamHeader --ignoreDummyR 40 --ignoreRepeatR 15 > clip.log &
    wait
    ls
    cd ..

}


function setup_clip_seed {
    SEED="F1"; MISMATCHES=1;
    if [ $LENGTH>=100 ]; then anchor=35
    elif [ $LENGTH==75 ]; then anchor=25
    else echo "UNSUPPORTED LEN"; fi 
}


function setup_clip_reads {
    CLIPDIR=$1; READS=$2
    mkdir -p $CLIPDIR 
    rm -f $CLIPDIR/gene_reads.txt; rm -f $CLIPDIR/genome_reads.txt
    for f in $READS; do 
        if [ ! -e $f ]; then echo "ERROR: FASTQ FILES NOT SUPPLIED" $f; exit; fi 
        rEXT=${f##*.} ; FNAME=$(readlink -m $f); FLINK=$CLIPDIR/$(basename $FNAME)
        if [ $rEXT != "fastq" ]; then echo "ERROR: FASTQ FILES NEEDED"; exit; fi 
        ln -fns $FNAME $FLINK; 
        printf $(basename $FNAME)"\n" >> $CLIPDIR/gene_reads.txt  
        printf $(basename $FNAME "."$rEXT)"_miss_GENES.fastq\n" >> $CLIPDIR/genome_reads.txt
    
    done
}

    

function setup_perm_reads {
    if [ $MAPDIR == "NONE" ]; then echo "ERROR: NEED AN OUTPUT DIRECTORY"; exit; fi 
   
     
    rm -f $MAPDIR/feature_reads.txt; rm -f $MAPDIR/genome_reads.txt; rm -f $MAPDIR/splice_reads.txt 
    for f in $READS; do
        if [ ! -e $f ]; then echo "ERROR: FASTQ FILES NOT SUPPLIED" $f; exit; fi 
        rEXT=${f##*.} ; FNAME=$(readlink -m $f); FLINK=$MAPDIR/$(basename $FNAME)
        if [ $rEXT != "fastq" ]; then echo "ERROR: FASTQ FILES NEEDED"; exit; fi 
        ln -fns $FNAME $FLINK; 
        printf $(basename $FNAME)"\n" >> $MAPDIR/feature_reads.txt  
        printf $(basename $FNAME "."$rEXT)"_miss_FEATURES.fastq\n" >> $MAPDIR/genome_reads.txt
        printf $(basename $FNAME "."$rEXT)"_miss_FEATURES_miss_GENOME.fastq\n" >> $MAPDIR/splice_reads.txt
    done
}

function setup_clip_refs {
    GENE_IDX="GENES_"$SEED"_"$LENGTH; GENOME_IDX="GENOME_"$SEED"_"$LENGTH

    if [ ! -z ${!GENE_IDX} ] && [ -f ${!GENE_IDX} ]; then GENES=${!GENE_IDX};
    elif [ -f $GENE_REF ]; then GENES=$GENE_REF
    else echo "INVALID REFS SUPPLIED"; exit; fi 
    if [ ! -z ${!GENOME_IDX} ]; then GENOME=${!GENOME_DIX}; else GENOME=$GENOME_REF; fi
    
    GENES_EXT=${GENES##*.}; GENOME_EXT=${GENOME##*.};
    ln -fns $GENES $CLIPDIR"/GENES."$GENES_EXT; myGENES="GENES."$GENES_EXT
    if [ ! -z $GENOME ]; then  ln -fns $GENOME $CLIPDIR/"GENOME."$GENOME_EXT; myGENOME="GENOME."$GENOME_EXT; fi
}



function setup_perm_refs {
    if [ $MAPDIR == "NONE" ]; then echo "ERROR: NEED AN OUTPUT DIRECTORY"; exit; fi 
        
    F_IDX="FEATURES_"$SEED"_"$LENGTH; G_IDX="GENOME_"$SEED"_"$LENGTH; S_IDX="SPLICE_"$SEED"_"$LENGTH
    
    if [ ! -z ${!F_IDX} ] && [ -f ${!F_IDX} ]; then FEATURES=${!F_IDX}; 
    elif [ -f $FEATURE_REF ]; then FEATURES=$FEATURE_REF
    else echo "INVALID REFERENCES SUPPLIED"; exit; fi 
    if [ ! -z ${!G_IDX} ]; then GENOME=${!G_IDX}; else GENOME=$GENOME_REF; fi 
    if [ ! -z ${!S_IDX} ]; then SPLICES=${!S_IDX}; else SPLICES=$SPLICE_REF; fi 

    FEATURE_EXT=${FEATURES##*.}; GENOME_EXT=${GENOME##*.}; SPLICE_EXT=${SPLICES##*.}  
    
    ln -fns $FEATURES $MAPDIR/"FEATURES."$FEATURE_EXT; myFEATURES="FEATURES."$FEATURE_EXT
    if [ ! -z $GENOME ]; then  ln -fns $GENOME $MAPDIR/"GENOME."$GENOME_EXT;  myGENOME="GENOME."$GENOME_EXT; fi 
    if [ ! -z $SPLICES ]; then ln -fns $SPLICES $MAPDIR/"SPLICES."$SPLICE_EXT;   mySPLICES="SPLICES."$SPLICE_EXT; fi 


}

function setup_perm_seeds {
    if [ $LENGTH -gt 64 ]; then  SEEDSCR=$(( ($MISMATCHES+1)/2 )); else SEEDSCR=$MISMATCHES; fi; SEED="F"$SEEDSCR
    
   }


    

function map_and_parse_reads_to_features {
    
    terminal_talk "Mapping reads to features..."
    valid_read_files=$(validate_reads feature_reads.txt val_reads.txt)  
    if [ -f FEATURES.log ] && [ $SKIP == "TRUE" ]; then terminal_talk "SKIPPING\n"
    elif [ $valid_read_files == 0 ]; then terminal_talk "NO READS\n"; 
    else perm $myFEATURES val_reads.txt --seed $SEED -v $MISMATCHES -B --printNM -u -s -T $LENGTH > FEATURES.log; check_progress $?; terminal_talk "Success\n"; fi 
    terminal_talk  "Parsing read alignments.."
    for i in *.mapping; do if [ -f $i".vis" ] && [ $SKIP == "TRUE" ]; then terminal_talk ".SKIP."; else parse_alignment.py $i --strandRule $STRANDRULE --tag $TAG > "$i".vis & fi done
    wait
    check_progress $?
    terminal_talk ".Sucess\n"
}




function map_and_parse_reads_to_genome {
    terminal_talk "Mapping reads to genome..."
    #validate_reads genome_reads.txt val_reads.txt  
    valid_read_files=$(validate_reads genome_reads.txt val_reads.txt)  
    if [ -f GENOME.log ] && [ $SKIP == "TRUE" ]; then terminal_talk "SKIPPING\n"
    elif [ $valid_read_files == 0 ]; then terminal_talk "NO READS\n"; 
    else perm $myGENOME val_reads.txt --seed $SEED -v $MISMATCHES -B --printNM -u -s -T $LENGTH --noSamHeader --outputFormat sam > GENOME.log; check_progress $?; terminal_talk "Success\n"; fi 
    
    terminal_talk "Parsing read alignments.."
    for i in *.sam; do if [ -f $i".vis" ] && [ $SKIP == "TRUE" ]; then terminal_talk ".SKIP."; else parse_alignment.py $i --tag $TAG  > "$i".vis & fi done
    wait
    check_progress $?
    terminal_talk ".Sucess\n"
}



function map_and_parse_reads_to_splices {
    terminal_talk "Mapping reads to splicesites..."; 
    valid_read_files=$(validate_reads splice_reads.txt val_reads.txt)  
    if [ -f SPLICE.log ] && [ $SKIP == "TRUE" ]; then terminal_talk "SKIPPING\n"
    elif [ $valid_read_files == 0 ]; then NOREADSZ=5; 
    else perm $mySPLICES val_reads.txt --seed $SEED -v $MISMATCHES -B --printNM -u -s -T $LENGTH > SPLICE.log; check_progress $?; terminal_talk "Success\nParsing read alignments.."; fi 
    
    for i in *.mapping; do  if [ -f $i".vis" ] && [ $SKIP == "TRUE" ]; then terminal_talk ".SKIP."; else parse_alignment.py $i --strandRule $STRANDRULE --tag $TAG > "$i".vis & fi done
    wait
    check_progress $?
    terminal_talk "..Sucess\n"
}

function validate_reads {
    # READS = 1 and outfile =2 
    rm -f $2
    passed=0
    for r in $(cat $1); do
        if [ -f $r ] && [ $(ls -l $(readlink -f $r) | awk '{print $5}') -gt 10 ]; then
           let passed+=1
           echo $r >> $2 
        fi
    done
    echo $passed
}




function terminal_talk {
    if [ $VOLUME != "SILENT" ]; then printf "$1"; fi }


function whisper {
    if [ $VOLUME == "SILENT" ]; then printf "$1"; fi }

function shout {
   if [ $VOLUME != "BAD" ]; then  printf "$1"; fi  }







   




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









