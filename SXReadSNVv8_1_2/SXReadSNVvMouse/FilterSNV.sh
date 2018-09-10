#!/bin/sh

APPDIR="/archive/project/oneclick/syhwah"

if [ "$M_SAMPLEID" = "" ];then SAMPLE_ID=$1;else SAMPLE_ID=$M_SAMPLEID;fi
if [ "$MIN_KSR" = "" ];then MIN_KSR=$2;else MIN_KSR=2;fi       #Min. Known Supporting Reads 
if [ "$3" != "" ];then MIN_NSR=$3;else MIN_NSR=2;fi            #Min. Novel Supporting Reads
if [ "$4" != "" ];then MAX_SR=$4;else MAX_SR=180;fi            #Max. Known/Novel Supporting Reads 
if [ "$5" != "" ];then MIN_RS=$5;else MIN_RS=15;fi             #Min. Read Strength
if [ "$6" != "" ];then MIN_RD=$6;else MIN_RD=4;fi              #Min. Read Density
if [ "$7" != "" ];then SNV_TYPE=$7;else SNV_TYPE=$M_TYPE; fi   #type=> S-Snp, I-Ins, D-Del


case "$SNV_TYPE" in
   "S")
       SNAME="snp"
       MIN_RS=20 
       ;;
   "I")
       SNAME="ins"
       ;;
   "D")
       SNAME="del"
       ;;
esac


CURR_DATE=$(date '+%y%m%d')
LIST=$SAMPLE_ID"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"
LIST_KNOWN=$SAMPLE_ID"_"$SNAME"_Known"
LIST_NOVEL=$SAMPLE_ID"_"$SNAME"_Novel"

#--------------------------------------------------------------------------------------------

PRG=$APPDIR"/SXParseAQSList/./SXParseAQSList"
echo $PRG $LIST $LIST_KNOWN".lst" $LIST_NOVEL".lst"
time $PRG $LIST $LIST_KNOWN".lst" $LIST_NOVEL".lst"

file_size=(`du -cb $LIST_KNOWN".lst"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $LIST_NOVEL".lst"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#--------------------------------------------------------------------------------------------

PRG=$APPDIR"/SXFilterSRB4Merge/./SXFilterSRB4Merge"
LIST_KNOWN_F1=$LIST_KNOWN"_F1.lst"

echo $PRG $LIST_KNOWN".lst" $LIST_KNOWN_F1 $MIN_KSR
time $PRG $LIST_KNOWN".lst" $LIST_KNOWN_F1 $MIN_KSR

file_size=(`du -cb $LIST_KNOWN_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_KNOWN_F2=$LIST_KNOWN"_F"$MAX_SR".lst"

echo $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $MAX_SR
time $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $MAX_SR 

file_size=(`du -cb $LIST_KNOWN_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

LIST_NOVEL_F1=$LIST_NOVEL"_F1.lst"

PRG=$APPDIR"/SXFilterSRRSRD/./SXFSRRSRD_SingleSample"
echo $PRG $LIST_NOVEL".lst" $MIN_NSR $MIN_RS $MIN_RD ">" $LIST_NOVEL_F1
eval time $PRG $LIST_NOVEL".lst" $MIN_NSR $MIN_RS $MIN_RD $LIST_NOVEL ">" $LIST_NOVEL_F1

file_size=(`du -cb $LIST_NOVEL_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_NOVEL_F2=$LIST_NOVEL"_F"$MAX_SR".lst"

echo $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $MAX_SR
time $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $MAX_SR

file_size=(`du -cb $LIST_NOVEL_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FILTER_LIST=$SAMPLE_ID"_"$SNAME"_Sample"

eval cat $LIST_KNOWN_F2 $LIST_NOVEL_F2 ">" $FILTER_LIST

file_size=(`du -cb $FILTER_LIST`)
if [ $file_size -eq 0 ]; then
   exit 1
fi



#------------------------------------------------------------------------------------------

M_FILTER_LIST=$FILTER_LIST #".lst"
#
#PRG=$APPDIR"/SXFilterSNPFalseHomo/./SXFilterSNPFalseHomo"
#
#echo $PRG $FILTER_LIST $M_FILTER_LIST "6.14"
#time $PRG $FILTER_LIST $M_FILTER_LIST "6.14"


