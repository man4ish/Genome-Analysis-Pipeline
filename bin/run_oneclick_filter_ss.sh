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
LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"

#--------------------------------------------------------------------------------------------

PRG="SXFilterSRB4Merge_SingleSample"
M_FILTER_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_"$SNAME"_Sample.SR"$M_MIN_SR

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $LIST $M_FILTER_LIST $M_MIN_SR >> $M_LOG_FILE 2>> $M_LOG_FILE

#------------------------------------------------------------------------------------------


