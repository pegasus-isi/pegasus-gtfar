#!/bin/bash



source /export/uec-gs1/knowles/analysis/tade/gtfar_source/optparse
export PATH=$PATH:/export/uec-gs1/knowles/analysis/tade/gtfar_source

#source /export/uec-gs1/knowles/analysis/jen_mapping/gtfar/gtfar_source/optparse
#export PATH=$PATH:/export/uec-gs1/knowles/analysis/jen_mapping/gtfar/gtfar_source


HG19=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19

#HGREF=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19/chr20.fa 
HGREF=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/hg19/hg19.txt

#DATA=/export/uec-gs1/knowles/analysis/tade/gtfar_runs/data

GTF=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/gencode16/gencode.v16.annotation.gtf
DATA=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/gencode16/annotation_data

INDEXES=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/indexes
GC16_IDXS=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/indexes/gc16
HG19_IDXS=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/indexes/gc16

#  INDEXES #

EXON_100_F1=$GC16_IDXS/50_F1_gc16_exonSeqs.index
EXON_100_F2=$GC16_IDXS/50_F2_gc16_exonSeqs.index
EXON_100_F3=$GC16_IDXS/50_F3_gc16_exonSeqs.index
EXON_100_F4=$GC16_IDXS/50_F4_gc16_exonSeqs.index

INTRON_100_F1=$GC16_IDXS/50_F1_gc16_intronSeqs.index
INTRON_100_F2=$GC16_IDXS/50_F2_gc16_intronSeqs.index
INTRON_100_F3=$GC16_IDXS/50_F3_gc16_intronSeqs.index
INTRON_100_F4=$GC16_IDXS/50_F4_gc16_intronSeqs.index

CATCLIP_40_F1=$GC16_IDXS/40_F1_gc16_catsOnly.index
GENECLIP_40_F1=$GC16_IDXS/40_F1_gc16_geneSeqs.index

HG19_100_F1=$HG19_IDXS/50_F1_hg19.index
HG19_100_F2=$HG19_IDXS/50_F2_hg19.index
HG19_100_F3=$HG19_IDXS/50_F3_hg19.index
HG19_100_F4=$HG19_IDXS/50_F4_hg19.index

HG19CLIP_40_F1=$HG19_IDXS/40_F1_hg19.index

KEY=$DATA/gc16.key
CATS=$DATA/gc16_catsOnly.fa
EXONS=$DATA/gc16_exonSeqs.fa
GENES=$DATA/gc16_geneSeqs.fa
INTRONS=$DATA/gc16_intronSeqs.fa
#GTF=$DATA/TEST19_20.gtf 


###  PROGRAMS:  annotate_references.py  call_mutations_from_locations.py  filter_fastq.py  modules  old  optparse  parse_sam_mapping.py



add_arg program 'annotation/filter-reads/iter-map' 

add_option prefix -p 'test' 'OutPut Prefix'
add_option reads -r 'reads.fq/fastq'   'read file in fq format'
add_option output -o 'output'   'output directory'
add_flag showcommands -h 1 'Print all executed commands to standard out' 


parseargs $@




READ_EXT=$(echo $READS | awk -F\. '{print $NF}')

if [ $READ_EXT != "fq" ] && [ $READ_EXT != "fastq" ]; then
    echo "A valid read file in fq format is required" 
fi

if [ -d $OUTPUT ]; then echo "WARNING: output directory ("$OUTPUT") already exists; Files may be overwritten"; 
else mkdir $OUTPUT; fi 




#####################################################  TOOLS            ###################################################################3

function smart_move {
    WC=(*$1)
    if [ -e ${WC[0]} ]; then
        mv *$1 $2/
    fi
}















######################################################  FEATURE MAPPING   #################################################################################

function first_feature_map {
       
    MYSEED=$1; MYSUBS=$2; MYREF=$3; MYREADS=$4; TYPE=$5; MYPREFIX=$6; MAPLOG=$7; SAY1=$8;
        
    if ! [ -z "$SAY1" ]; then echo  "Mapping reads to "$SAY1" (Iteration 1 )..."; fi 
    perm $MYREF $MYREADS --seed $MYSEED -v $MYSUBS -B -o $MYPREFIX.mapping --printNM -u "$MYPREFIX"_miss.fastq -s "$MYPREFIX".index >> $MAPLOG
    if ! [ -z "$SAY1" ]; then echo "Parsing alignment..."; fi 
    parse_alignment.py $TYPE $MYPREFIX.mapping -k $KEY -p "$MYPREFIX"_parse 
    if ! [ -z "$SAY1" ]; then echo "Determining mutations"...; fi 
    call_mutations.py "$MYPREFIX"_parse_gene.srt $GTF $TYPE -g $HG19 -p "$MYPREFIX"_iter 
    if ! [ -z "$SAY1" ]; then echo "Complete"; echo ""; fi 
    }

function second_feature_map {   
    MYSEED=$1; MYSUBS=$2; MYREF=$3; MYREADS=$4; TYPE=$5; MYPREFIX=$6; MAPLOG=$7; SAY1=$8
    if ! [ -z "$SAY1" ]; then echo "Mapping reads to "$SAY1" (Iteration 2 )..."; fi 
    perm $MYREF $MYREADS --seed $MYSEED -v $MYSUBS -B -o $MYPREFIX.mapping --printNM -u "$MYPREFIX"_miss.fastq -s "$MYPREFIX".index >> $MAPLOG
    if ! [ -z "$SAY1" ]; then echo "Parsing alignment..."; fi 
    parse_alignment.py $TYPE $MYPREFIX.mapping -k $KEY -p "$MYPREFIX"_parse 
    }

function iterative_feature_analysis {

    SEED1=$1; SUB1=$2; SEED2=$3; SUB2=$4; MAPLOG=$5; SAY=$6
    # ROUND 1 # 
    first_feature_map $SEED1 $SUB1 $EXON_100_F1 $READS --exonic $OUTPUT/exonic1 $MAPLOG 'exonic seqs'
    first_feature_map $SEED1 $SUB1 $INTRON_100_F1 $OUTPUT/exonic1_miss.fastq --intronic $OUTPUT/intronic1 $MYLOG "intronic seqs"
    # ROUND 2 # 
    second_feature_map $SEED2 $SUB2 $OUTPUT/exonic1_iter_exonSeqs.fa $OUTPUT/intronic1_miss.fastq --exonic $OUTPUT/exonic2 $MYLOG 'exonic seqs' 
    second_feature_map $SEED2 $SUB2 $OUTPUT/intronic1_iter_intronSeqs.fa $OUTPUT/exonic2_miss.fastq --intronic $OUTPUT/intronic2 $MYLOG 'intron seqs' 
    # COMBINE ANALYSIS # 

    if ! [ -z "$SAY" ]; then echo "Combining/mergeing iterative results..."; fi 
    sort -m -k14n,14 -k6,6 -k8n,8 $OUTPUT/exonic*.srt >   $OUTPUT/exonic_iter_gene.srt 
    sort -m -k14n,14 -k6,6 -k8n,8 $OUTPUT/intronic*.srt > $OUTPUT/intronic_iter_gene.srt  
    if ! [ -z "$SAY" ]; then echo "Determining iterative mutations..."; fi 
    call_mutations.py $OUTPUT/exonic_iter_gene.srt $GTF --exonic -g $HG19 -p $OUTPUT/exonic_iterative 
    call_mutations.py $OUTPUT/intronic_iter_gene.srt $GTF --intronic -g $HG19 -p $OUTPUT/intronic_iterative 

    }

###########################################################################################################################################################################




function double_genome_map { 
    SEED=$1; SUBS=$2; MYREADS=$3; MYPREFIX=$4; MAPLOG=$5; SAY=$6; 
    if ! [ -z "$SAY" ]; then echo "Mapping reads to human genome..."; fi 
    perm $HG19_100_F2 $MYREADS --seed $SEED -v $SUBS -B -o $MYPREFIX.mapping --printNM -u "$MYPREFIX"_miss.fastq >> $MAPLOG  
    perm $HG19_100_F2 $MYREADS --seed $SEED -v $SUBS -B -o $MYPREFIX.sam > "$MYPREFIX"_sam.log
    }



function multi_gapped_alignment {

    REF1=$1; REF2=$2; REF3=$3; MYREADS=$4; MYPREFIX=$5; MAPLOG=$6; SAY=$7
    
    if ! [ -z "$SAY" ]; then echo "Gapped mapping to concatenated exons"; fi 
    clipR $REF1 $MYREADS --anchorL 40 --seed F1 -v 1 --ignoreRepeatR 10 --ignoreDummyR 40 -u "$MYPREFIX"_cat_miss.fastq -o "$MYPREFIX"_cat.sam --noSamHeader >> $MAPLOG  
    if ! [ -z "$SAY" ]; then echo "Parsing gapped exonic alignment"; fi 
    parse_alignment.py "$MYPREFIX"_cat.sam --gapped -k $KEY -p "$MYPREFIX"_cat_parse    
    if ! [ -z "$SAY" ]; then echo "Gapped mapping to gene sequences"; fi 
    clipR $REF2 "$MYPREFIX"_cat_miss.fastq --anchorL 40 --seed F1 -v 1 --ignoreRepeatR 10 --ignoreDummyR 40 -u "$MYPREFIX"_gene_miss.fastq -o "$MYPREFIX"_gene.sam --noSamHeader >> $MAPLOG 
    if ! [ -z "$SAY" ]; then echo "Parsing gapped gene alignment"; fi 
    parse_alignment.py "$MYPREFIX"_gene.sam --gapped -k $KEY -p "$MYPREFIX"_gene_parse    
    if ! [ -z "$SAY" ]; then echo "Gapped mapping to human genome"; fi 
    clipR $REF3 "$MYPREFIX"_gene_miss.fastq --anchorL 40 --seed F1 -v 1 --ignoreRepeatR 10 --ignoreDummyR 40 -u "$MYPREFIX"_hg19_miss.fastq -o "$MYPREFIX"_hg19.sam --noSamHeader >> $MAPLOG  
    if ! [ -z "$SAY" ]; then echo "Parsing gapped genomic alignment"; fi 
    parse_alignment.py "$MYPREFIX"_hg19.sam --hg19 -k $KEY -p "$MYPREFIX"_hg19_parse    
    }




if [ $PROGRAM == "iter-map" ] || [ $PROGRAM == "ITER-MAP" ]; then 

        # ---- MAP STEP -----  #

        MYLOG=$OUTPUT/MAPPING.stats 

        # ITERATIVE FEATURE ALIGNMENT # 

        echo "Pipeline Begins Sample: " $READS


        iterative_feature_analysis F1 2 F3 5 $MYLOG GO

        # GENOMIC ALIGNMENT (TWO FILETYPES) #
        double_genome_map F2 4 $OUTPUT/intronic2_miss.fastq $OUTPUT/hg19 $MYLOG GO 
        # GAPPED ALIGNMENT - MULTIPLE WAYS #

        multi_gapped_alignment $OUTPUT/exonic_iterative_catsOnly.fa $OUTPUT/exonic_iterative_geneSeqs.fa $HG19CLIP_40_F1 $OUTPUT/hg19_miss.fastq $OUTPUT/GAPS $MYLOG GO 
        cp $OUTPUT/GAPS_hg19_miss.fastq $OUTPUT/"$PREFIX_"remainingReads.fq
        echo "MAPPING COMPLETE"
        echo ""

        parse_log.py $MYLOG 

        if  [ ! -d $OUTPUT/results ]; then mkdir $OUTPUT/results;  fi
        if  [ ! -d $OUTPUT/tmp     ]; then mkdir $OUTPUT/tmp;  fi

        ## CLEANUP ROUND ##

        cd $OUTPUT 
        smart_move .mapping tmp
        smart_move .fa      tmp
        smart_move .log     tmp
        smart_move .loc     tmp
        smart_move .fastq  tmp
        cd ..
        # MERGE EXPRESSION, ETC, ETC, ETC # --- ALSO YOU COULD PARSE THE GENOME BETTER  
        
        cat $OUTPUT/*gene.cnts | sort -k1,1 -k5,5 | awk '{if ($5!=L) {if (NR!=1) print A,B,C,"|",L,X,Y,Z; A=$1;B=$2;C=$3;L=$5;X=$6;Y=$7;Z=$8} else {X+=$6;Y+=$7;Z+=$8}}' > $OUTPUT/results/$PREFIX.gene.cnts
    fi



exit




        ###  MERGING FILES FOR MUTATION CALLING ###

   #     sort -m -k14n,14 -k6,6 -k8n,8 $OUTPUT/ex*.srt > $OUTPUT/exJoin_parse_gene.srt 
   #     sort -m -k14n,14 -k6,6 -k8n,8 $OUTPUT/int*.srt > $OUTPUT/intJoin_parse_gene.srt 
   #     call_mutations.py output/exJoin_parse_gene.srt  $GTF --exonic    -g $HG19 -p $OUTPUT/exonic_iterative
   #     call_mutations.py output/intJoin_parse_gene.srt $GTF --intronic -g $HG19 -p $OUTPUT/intronic_iterative

        ### USED MERGE FILES FOR GAPPED ALIGNMENT ###

        echo "Begin Gapped Alignment"

   #     clipR $OUTPUT/exonic_iterative_catsOnly.fa $OUTPUT/miss_round2.fq -o $OUTPUT/gapped_exonic.sam --noSamHeader -e --anchorL 22 --seed F1 -e -v 4 -s $OUTPUT/clipSeq.idx -u $OUPUT/clipMiss.fq --ignoreRepeatR 15 --ignoreDummyR 40 > $OUTPUT/clip.log


#ex1_iter_catsOnly.fa  ex1_iter_exonSeqs.fa  int1_iter_intronSeqs.fa



    fi
exit 





if [ -d reads ]; then echo "WARNING: reads directory already exists; Files may be overwritten"
else mkdir reads; fi
if [ -d refs ]; then echo "WARNING: refs directory  already exists; Files may be overwritten"
else mkdir refs; fi 
if [ -d maps ]; then echo "WARNING: maps directory  already exists; Files may be overwritten"
else mkdir maps; fi
if [ -d logs ]; then echo "WARNING: maps directory  already exists; Files may be overwritten"
else mkdir logs; fi


