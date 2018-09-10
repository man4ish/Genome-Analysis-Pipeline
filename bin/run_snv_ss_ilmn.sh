#!/bin/sh

if [ -z "$SXDB" ] ; then
   SXDB="sprt_prtn_upkb-57.14_sb3.0.9"
fi

#checked the availability of environment variables for synasearch
. "check_synasearch_env.sh"

if [ $? -ne "0" ]; then
   exit 1
fi

ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <Project Home Dir> <Sample-ID> <Genome Ver -eg. 37.1> <Coverage> \
      <Human Ref. files path> [synaworks_resources] [QScore Cutoff]" 

if [ $# -lt 5 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi



#------------------------------------------------------------------------------------------------
export M_PROJECTDIR=$1
export M_SAMPLEID=$2
export M_GENOME=$3
export M_MEAN=$4
export M_HREF_FILEPATH=$5
export M_READSDIR=$M_PROJECTDIR"/density"
export M_OUTPUT_DIR=$M_PROJECTDIR"/snv/tmp"
export M_FINALOUT_DIR=$M_PROJECTDIR"/snv"
export M_PROJ_SAMPLEID=$M_SAMPLEID

CURR_DATE=$(date '+%y%m%d')
export M_LOG_FILE=$M_OUTPUT_DIR"/snv_"$M_SAMPLEID"_"$CURR_DATE".log"

if [ $# -gt 5 ]; then
   export M_REFDIR=$6
else
   export M_REFDIR="/data/dbref/synaworks/resources/genome_37.1" 
fi
FDBSNP=$M_REFDIR"/dbSNP_131"
FDBINS=$M_REFDIR"/dbINS_131"
FDBDEL=$M_REFDIR/"dbDEL_131" 

if [ $# -gt 6 ]; then
   QSCUTOFF=$7
else
   QSCUTOFF=64 
fi

mkdir -p $M_OUTPUT_DIR

if [ $M_GENOME = "37.1" ]; then
   PRG_OPTION="-gLa"
else
   PRG_OPTION="-gLaO" 
fi

#----------------------------------------------------------------------------------------------------------------------------

PRG="SXReadSNV_Ilm"
RAW_FILE=$M_PROJECTDIR"/snv/"$M_SAMPLEID".snv"
RAW_FILE_FILTERED=$RAW_FILE
#".filtered"

eval cat $M_PROJECTDIR"/snv/"$M_SAMPLEID".snv_*" > $RAW_FILE

echo "" >> $M_LOG_FILE
#eval /usr/bin/time --verbose $PRG -f $RAW_FILE -fS $QSCUTOFF 3 100 $RAW_FILE_FILTERED >> $M_LOG_FILE 2>> $M_LOG_FILE

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG -f $RAW_FILE_FILTERED $PRG_OPTION $FDBSNP $FDBINS $FDBDEL $M_READSDIR 1 10000 0 0 $M_SAMPLEID \
                             $M_GENOME $QSCUTOFF $M_HREF_FILEPATH $M_OUTPUT_DIR $M_PROJ_SAMPLEID -c >> $M_LOG_FILE 2>> $M_LOG_FILE
#----------------------------------------------------------------------------------------------------------------------------

CONFIG=$M_PROJECTDIR"/conf/snv.conf"

SNP_MIN_SR=3
SNP_MIN_CSCORE=0.6

INS_MIN_SR=3
INS_MIN_CSCORE=0.7

DEL_MIN_SR=3
DEL_MIN_CSCORE=0.7


if [ -f $CONFIG ]; then
   file_size=(`du -cb $CONFIG`)
   if [ $file_size -ne 0 ]; then
      tmp=$(grep SNP_MIN_SR $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')    
      if [ -n "$tmp" ]; then SNP_MIN_SR=$tmp; fi

      tmp=$(grep SNP_MIN_CSCORE $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')
      if [ -n "$tmp" ]; then SNP_MIN_CSCORE=$tmp; fi

      tmp=$(grep INS_MIN_SR $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')
      if [ -n "$tmp" ]; then INS_MIN_SR=$tmp; fi

      tmp=$(grep INS_MIN_CSCORE $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')
      if [ -n "$tmp" ]; then INS_MIN_CSCORE=$tmp; fi

      tmp=$(grep DEL_MIN_SR $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')
      if [ -n "$tmp" ]; then DEL_MIN_SR=$tmp; fi

      tmp=$(grep DEL_MIN_CSCORE $CONFIG | cut -d= -f2 | sed 's/^ *//g;s/ *$//g')
      if [ -n "$tmp" ]; then DEL_MIN_CSCORE=$tmp; fi
   fi
fi


export M_MIN_SR=$SNP_MIN_SR #3
export M_MIN_CSCORE=$SNP_MIN_CSCORE #0.7

echo "" >> $M_LOG_FILE
echo "SNP_MIN_SR="$SNP_MIN_SR >> $M_LOG_FILE 
echo "SNP_MIN_CSCORE="$SNP_MIN_CSCORE >> $M_LOG_FILE 
echo "INS_MIN_SR="$INS_MIN_SR >> $M_LOG_FILE
echo "INS_MIN_CSCORE="$INS_MIN_CSCORE >> $M_LOG_FILE
echo "DEL_MIN_SR="$DEL_MIN_SR >> $M_LOG_FILE
echo "DEL_MIN_CSCORE="$DEL_MIN_CSCORE >> $M_LOG_FILE
echo "" >> $M_LOG_FILE

export M_FILTER_LIST=""     #value will be defined in ProcSnp.sh
export M_MIN_PSR="1"        #min supporting reads to be printed in final stats report
export M_MAX_PSR="500"      #max supporting reads to be printed in final stats report
export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del  


#Process SNP List ...
. "run_oneclick_filter_ss.sh"

#-----------------------------------------------------------------------------------------------

M_TYPE="S"

. "run_oneclick_annotate_ss_snp_ilmn.sh"

#-----------------------------------------------------------------------------------------------

#Process INS List ...

M_TYPE="I"
M_SAMPLENO=1
M_MIN_SR=$INS_MIN_SR
M_MIN_CSCORE=$INS_MIN_CSCORE

. "run_oneclick_filter_ss.sh"
. "run_oneclick_annotate_ss_ilmn.sh"

#Process DEL List ...

M_TYPE="D"
M_MIN_SR=$DEL_MIN_SR 
M_MIN_CSCORE=$DEL_MIN_CSCORE

. "run_oneclick_filter_ss.sh"
. "run_oneclick_annotate_ss_ilmn.sh"
