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

#--------------------------------------------------------------------------------------------

if [ "$M_TYPE" != "S" ]; then
   PRG=$M_APPDIR"/SXMergePoly/./SXMergePoly"
   echo $PRG $LIST $M_TYPE $LIST".merge"
   time $PRG $LIST $M_TYPE $LIST".merge" 

   #file_size=(`du -cb $LIST".merge")
   #if [ $file_size -eq 0 ]; then
   #   exit 1
   #fi

   mv $LIST $LIST".org"
   mv $LIST".merge" $LIST
fi


#--------------------------------------------------------------------------------------------

LIST_KNOWN=$M_SAMPLEID"_"$SNAME"_Known"
LIST_NOVEL=$M_SAMPLEID"_"$SNAME"_Novel"

#--------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXParseAQSList/./SXParseAQSList"
echo $PRG $LIST $LIST_KNOWN $LIST_NOVEL
time $PRG $LIST $LIST_KNOWN $LIST_NOVEL


PRG=$M_APPDIR"/SXFilterSRB4Merge/./SXFilterSRB4Merge_SingleSample"
LIST_KNOWN_F1=$LIST_KNOWN"_F1"

echo $PRG $LIST_KNOWN $LIST_KNOWN_F1 $M_MIN_KSR
time $PRG $LIST_KNOWN $LIST_KNOWN_F1 $M_MIN_KSR

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_KNOWN_F2=$LIST_KNOWN"_F"$M_MAX_SR

echo $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR
time $PRG $LIST_KNOWN_F1 $LIST_KNOWN_F2 $M_MAX_SR 

#------------------------------------------------------------------------------------------
#This implementation is only meant for NA19240
#

if [ "$M_TYPE" = "S" ]; then
   PRG=$M_APPDIR"/SXFilterSR/./SXFilterSR_Grep_SingleSample"
   LIST_KNOWN_SR1_RS3=$LIST_KNOWN"_SR1_RS3"

   echo $PRG $LIST_KNOWN 1 3 $LIST_KNOWN_SR1_RS3
   time $PRG $LIST_KNOWN 1 3 $LIST_KNOWN_SR1_RS3

fi

#------------------------------------------------------------------------------------------

LIST_NOVEL_F1=$LIST_NOVEL"_F1"

PRG=$M_APPDIR"/SXFilterSRRSRD/./SXFSRRSRD_SingleSample"
echo $PRG $LIST_NOVEL $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1
eval time $PRG $LIST_NOVEL $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD ">" $LIST_NOVEL_F1

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT_SingleSample"
LIST_NOVEL_F2=$LIST_NOVEL"_F"$M_MAX_SR

echo $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR
time $PRG $LIST_NOVEL_F1 $LIST_NOVEL_F2 $M_MAX_SR

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

#------------------------------------------------------------------------------------------
