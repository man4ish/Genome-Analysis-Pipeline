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
echo $PRG $LIST $LIST_KNOWN $LIST_NOVEL
time $PRG $LIST $LIST_KNOWN $LIST_NOVEL

file_size=(`du -cb $LIST_KNOWN`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $LIST_NOVEL`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#--------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterSRB4Merge/./SXFilterSRB4Merge_SingleSample"
LIST_KNOWN_F1=$LIST_KNOWN"_F1"

echo $PRG $LIST_KNOWN $LIST_KNOWN_F1 $M_MIN_KSR
time $PRG $LIST_KNOWN $LIST_KNOWN_F1 $M_MIN_KSR

file_size=(`du -cb $LIST_KNOWN_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_KNOWN_F2=$LIST_KNOWN"_F"$M_MAX_SR

echo $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR
time $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR 

file_size=(`du -cb $LIST_KNOWN_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
#This implementation is only meant for NA19240
#

if [ "$M_TYPE" = "S" ]; then
   PRG=$M_APPDIR"/SXFilterSR/./SXFilterSR_Grep_SingleSample"
   LIST_KNOWN_SR1_RS3=$LIST_KNOWN"_SR1_RS3"

   echo $PRG $LIST_KNOWN 1 3 $LIST_KNOWN_SR1_RS3
   time $PRG $LIST_KNOWN 1 3 $LIST_KNOWN_SR1_RS3

   file_size=(`du -cb $LIST_KNOWN_SR1_RS3`)
   if [ $file_size -eq 0 ]; then
      exit 1
   fi
fi

#------------------------------------------------------------------------------------------

LIST_NOVEL_F1=$LIST_NOVEL"_F1"

PRG=$M_APPDIR"/SXFilterSRRSRD/./SXFSRRSRD_SingleSample"
echo $PRG $LIST_NOVEL $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1
eval time $PRG $LIST_NOVEL $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1

file_size=(`du -cb $LIST_NOVEL_F1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_NOVEL_F2=$LIST_NOVEL"_F"$M_MAX_SR

echo $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR
time $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR

file_size=(`du -cb $LIST_NOVEL_F2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

M_FILTER_LIST=$M_SAMPLEID"_"$SNAME"_Sample"


if [ "$M_TYPE" = "S" ]; then
   eval cat $LIST_KNOWN_F2 $LIST_KNOWN_SR1_RS3 $LIST_NOVEL_F2 ">" $M_FILTER_LIST

   PRG=$M_APPDIR"/SXCalculateRSGraph2d/./SXCalculateRSGraph2d_SingleSample"
   echo $PRG $FILTER_LIST $M_FILTER_LIST".Graph"
   time $PRG $FILTER_LIST $M_FILTER_LIST".Graph"  
else
   eval cat $LIST_KNOWN_F2 $LIST_NOVEL_F2 ">" $M_FILTER_LIST
fi

file_size=(`du -cb $M_FILTER_LIST`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
