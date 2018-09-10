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
LIST_SAMPLE=$M_SAMPLEID"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"

PRG=$M_APPDIR"/SXMergePoly/./SXMergePoly"

file_size=(`du -cb $LIST_SAMPLE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo $PRG $LIST_SAMPLE $M_TYPE $LIST_SAMPLE".merge"
time $PRG $LIST_SAMPLE $M_TYPE $LIST_SAMPLE".merge"

file_size=(`du -cb $LIST_SAMPLE".merge"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

mv $LIST_SAMPLE $LIST_SAMPLE".org"
