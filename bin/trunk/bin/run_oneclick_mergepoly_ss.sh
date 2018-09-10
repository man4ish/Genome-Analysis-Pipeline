#!/bin/sh

case "$M_TYPE" in
   "I")
       SNAME="ins"
       ;;
   "D")
       SNAME="del"
       ;;
esac

CURR_DATE=$(date '+%y%m%d')
LIST_SAMPLE=$M_OUTPUT_DIR"/"$M_SAMPLEID"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"

PRG="SXMergePoly"

file_size=(`du -cb $LIST_SAMPLE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $LIST_SAMPLE $M_TYPE $LIST_SAMPLE".merge" >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $LIST_SAMPLE".merge"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

eval /usr/bin/time --verbose mv $LIST_SAMPLE $LIST_SAMPLE".org" 2>> $M_LOG_FILE
