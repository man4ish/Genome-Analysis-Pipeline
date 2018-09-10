#!/bin/sh

case "$M_TYPE" in
   "S")
       SNAME="snp"
       ;;
   "I")
       SNAME="ins"
       ;;
   "D")
       SNAME="del"
       ;;
esac


CURR_DATE=$(date '+%y%m%d')
LIST=$M_SAMPLEID"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"
LIST_KNOWN=$M_SAMPLEID"_"$SNAME"_Known"
LIST_NOVEL=$M_SAMPLEID"_"$SNAME"_Novel"

#--------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXParseAQSList/./SXParseAQSList"
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

PRG=$M_APPDIR"/SXFilterSRB4Merge/./SXFilterSRB4Merge"
LIST_KNOWN_F1=$LIST_KNOWN"_F1.lst"

echo $PRG $LIST_KNOWN".lst" $LIST_KNOWN_F1 $M_MIN_KSR
time $PRG $LIST_KNOWN".lst" $LIST_KNOWN_F1 $M_MIN_KSR

file_size=(`du -cb $LIST_KNOWN_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_KNOWN_F2=$LIST_KNOWN"_F"$M_MAX_SR".lst"

echo $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR
time $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR 

file_size=(`du -cb $LIST_KNOWN_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

LIST_NOVEL_F1=$LIST_NOVEL"_F1.lst"

PRG=$M_APPDIR"/SXFilterSRRSRD/./SXFSRRSRD_SingleSample"
echo $PRG $LIST_NOVEL".lst" $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1
eval time $PRG $LIST_NOVEL".lst" $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1

file_size=(`du -cb $LIST_NOVEL_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_NOVEL_F2=$LIST_NOVEL"_F"$M_MAX_SR".lst"

echo $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR
time $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR

file_size=(`du -cb $LIST_NOVEL_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FILTER_LIST=$M_SAMPLEID"_"$SNAME"_Sample"

eval cat $LIST_KNOWN_F2 $LIST_NOVEL_F2 ">" $FILTER_LIST

file_size=(`du -cb $FILTER_LIST`)
if [ $file_size -eq 0 ]; then
   exit 1
fi



#------------------------------------------------------------------------------------------

M_FILTER_LIST=$FILTER_LIST #".lst"

PRG=$M_APPDIR"/SXCalculateRSGraph2d/./SXCalculateRSGraph2d_SingleSample"
echo $PRG $FILTER_LIST $M_FILTER_LIST".Graph"
time $PRG $FILTER_LIST $M_FILTER_LIST".Graph"

#
#PRG=$APPDIR"/SXFilterSNPFalseHomo/./SXFilterSNPFalseHomo"
#
#echo $PRG $FILTER_LIST $M_FILTER_LIST "6.14"
#time $PRG $FILTER_LIST $M_FILTER_LIST "6.14"


