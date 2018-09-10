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
LIST_SAMPLE_1=$M_SAMPLEID_1"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"
LIST_KNOWN_SAMPLE_1=$M_SAMPLEID_1"_"$SNAME"_Known"
LIST_NOVEL_SAMPLE_1=$M_SAMPLEID_1"_"$SNAME"_Novel"

#-------------------------------------------------------------------------------------------
#Merge INS and DEL polybase Sample1

#PRG=$M_APPDIR"/SXMergePoly/./SXMergePoly"

#file_size=(`du -cb $LIST_SAMPLE_1`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#echo $PRG $LIST_SAMPLE_1 $M_TYPE $LIST_SAMPLE_1".merge"
#time $PRG $LIST_SAMPLE_1 $M_TYPE $LIST_SAMPLE_1".merge"

#file_size=(`du -cb $LIST_SAMPLE_1".merge"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#mv $LIST_SAMPLE_1 $LIST_SAMPLE_1".org"
#mv $LIST_SAMPLE_1".merge" $LIST_SAMPLE_1


#------------------------------------------------------------------------------------------------
PRG=$M_APPDIR"/SXParseAQSList/./SXParseAQSList"

echo ""
echo $PRG $LIST_SAMPLE_1 $LIST_KNOWN_SAMPLE_1 $LIST_NOVEL_SAMPLE_1
time $PRG $LIST_SAMPLE_1 $LIST_KNOWN_SAMPLE_1 $LIST_NOVEL_SAMPLE_1

file_size=(`du -cb $LIST_KNOWN_SAMPLE_1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $LIST_NOVEL_SAMPLE_1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#--------------------------------------------------------------------------------------------

LIST_SAMPLE_2=$M_SAMPLEID_2"_"$SNAME"_c1_aqs_"$CURR_DATE".lst"
LIST_KNOWN_SAMPLE_2=$M_SAMPLEID_2"_"$SNAME"_Known"
LIST_NOVEL_SAMPLE_2=$M_SAMPLEID_2"_"$SNAME"_Novel"

#Merge INS and DEL polybase Sample2

#PRG=$M_APPDIR"/SXMergePoly/./SXMergePoly"

#file_size=(`du -cb $LIST_SAMPLE_2`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#echo $PRG $LIST_SAMPLE_2 $M_TYPE $LIST_SAMPLE_2".merge"
#time $PRG $LIST_SAMPLE_2 $M_TYPE $LIST_SAMPLE_2".merge"

#file_size=(`du -cb $LIST_SAMPLE_2".merge"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#mv $LIST_SAMPLE_2 $LIST_SAMPLE_2".org"
#mv $LIST_SAMPLE_2".merge" $LIST_SAMPLE_2

#------------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXParseAQSList/./SXParseAQSList"

echo ""
echo $PRG $LIST_SAMPLE_2 $LIST_KNOWN_SAMPLE_2 $LIST_NOVEL_SAMPLE_2
time $PRG $LIST_SAMPLE_2 $LIST_KNOWN_SAMPLE_2 $LIST_NOVEL_SAMPLE_2

file_size=(`du -cb $LIST_KNOWN_SAMPLE_2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $LIST_NOVEL_SAMPLE_2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#--------------------------------------------------------------------------------------------

SAMPLE1_SAMPLE2_KNOWN_COMMON=$M_SAMPLEID_1"_"$M_SAMPLEID_2"_"$SNAME"_Known_Common"
SAMPLE1_KNOWN_UNIQUE=$M_SAMPLEID_1"_"$SNAME"_Known_Unique"
SAMPLE2_KNOWN_UNIQUE=$M_SAMPLEID_2"_"$SNAME"_Known_Unique"

PRG=$M_APPDIR"/SXMerge/./SXMerge"

echo ""
echo $PRG $LIST_KNOWN_SAMPLE_1 $LIST_KNOWN_SAMPLE_2 $SAMPLE1_SAMPLE2_KNOWN_COMMON $SAMPLE1_KNOWN_UNIQUE $SAMPLE2_KNOWN_UNIQUE
time $PRG $LIST_KNOWN_SAMPLE_1 $LIST_KNOWN_SAMPLE_2 $SAMPLE1_SAMPLE2_KNOWN_COMMON $SAMPLE1_KNOWN_UNIQUE $SAMPLE2_KNOWN_UNIQUE

file_size=(`du -cb $SAMPLE1_SAMPLE2_KNOWN_COMMON`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $SAMPLE1_KNOWN_UNIQUE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $SAMPLE2_KNOWN_UNIQUE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#--------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterSR/./SXFilterSR"

echo ""
echo $PRG $SAMPLE1_SAMPLE2_KNOWN_COMMON $SAMPLE1_SAMPLE2_KNOWN_COMMON"_f1" $M_MIN_COMMON_KSR
time $PRG $SAMPLE1_SAMPLE2_KNOWN_COMMON $SAMPLE1_SAMPLE2_KNOWN_COMMON"_f1" $M_MIN_COMMON_KSR

file_size=(`du -cb $SAMPLE1_SAMPLE2_KNOWN_COMMON"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_KNOWN_UNIQUE $SAMPLE1_KNOWN_UNIQUE"_f1" $M_MIN_KSR
time $PRG $SAMPLE1_KNOWN_UNIQUE $SAMPLE1_KNOWN_UNIQUE"_f1" $M_MIN_KSR

file_size=(`du -cb $SAMPLE1_KNOWN_UNIQUE"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_KNOWN_UNIQUE $SAMPLE2_KNOWN_UNIQUE"_f1" $M_MIN_KSR2
time $PRG $SAMPLE2_KNOWN_UNIQUE $SAMPLE2_KNOWN_UNIQUE"_f1" $M_MIN_KSR2

file_size=(`du -cb $SAMPLE2_KNOWN_UNIQUE"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT"

SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED=$SAMPLE1_SAMPLE2_KNOWN_COMMON"_f"$M_MAX_SR
SAMPLE1_KNOWN_UNIQUE_FILTERED=$SAMPLE1_KNOWN_UNIQUE"_f"$M_MAX_SR
SAMPLE2_KNOWN_UNIQUE_FILTERED=$SAMPLE2_KNOWN_UNIQUE"_f"$M_MAX_SR

echo ""
echo $PRG $SAMPLE1_SAMPLE2_KNOWN_COMMON"_f1" $SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED $M_MAX_SR
time $PRG $SAMPLE1_SAMPLE2_KNOWN_COMMON"_f1" $SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED $M_MAX_SR 

file_size=(`du -cb $SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_KNOWN_UNIQUE"_f1" $SAMPLE1_KNOWN_UNIQUE_FILTERED $M_MAX_SR
time $PRG $SAMPLE1_KNOWN_UNIQUE"_f1" $SAMPLE1_KNOWN_UNIQUE_FILTERED $M_MAX_SR

file_size=(`du -cb $SAMPLE1_KNOWN_UNIQUE_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_KNOWN_UNIQUE"_f1" $SAMPLE2_KNOWN_UNIQUE_FILTERED $M_MAX_SR
time $PRG $SAMPLE2_KNOWN_UNIQUE"_f1" $SAMPLE2_KNOWN_UNIQUE_FILTERED $M_MAX_SR

file_size=(`du -cb $SAMPLE2_KNOWN_UNIQUE_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

SAMPLE1_SAMPLE2_NOVEL_COMMON=$M_SAMPLEID_1"_"$M_SAMPLEID_2"_"$SNAME"_Novel_Common"
SAMPLE1_NOVEL_UNIQUE=$M_SAMPLEID_1"_"$SNAME"_Novel_Unique"
SAMPLE2_NOVEL_UNIQUE=$M_SAMPLEID_2"_"$SNAME"_Novel_Unique"

PRG=$M_APPDIR"/SXMerge/./SXMerge"

echo ""
echo $PRG $LIST_NOVEL_SAMPLE_1 $LIST_NOVEL_SAMPLE_2 $SAMPLE1_SAMPLE2_NOVEL_COMMON $SAMPLE1_NOVEL_UNIQUE $SAMPLE2_NOVEL_UNIQUE
time $PRG $LIST_NOVEL_SAMPLE_1 $LIST_NOVEL_SAMPLE_2 $SAMPLE1_SAMPLE2_NOVEL_COMMON $SAMPLE1_NOVEL_UNIQUE $SAMPLE2_NOVEL_UNIQUE

file_size=(`du -cb $SAMPLE1_SAMPLE2_NOVEL_COMMON`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $SAMPLE1_KNOWN_UNIQUE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

file_size=(`du -cb $SAMPLE2_KNOWN_UNIQUE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterSRRSRD/./SXFSRRSRD"

echo ""
echo $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON $M_MIN_COMMON_NSR $M_MIN_NRS $M_MIN_COMMON_NRD $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f1"
time $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON $M_MIN_COMMON_NSR $M_MIN_NRS $M_MIN_COMMON_NRD $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f1"

file_size=(`du -cb $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_NOVEL_UNIQUE $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD $SAMPLE1_NOVEL_UNIQUE"_f1"
time $PRG $SAMPLE1_NOVEL_UNIQUE $M_MIN_NSR $M_MIN_NRS $M_MIN_NRD $SAMPLE1_NOVEL_UNIQUE"_f1"

file_size=(`du -cb $SAMPLE1_NOVEL_UNIQUE"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_NOVEL_UNIQUE $M_MIN_NSR2 $M_MIN_NRS $M_MIN_NRD2 $SAMPLE2_NOVEL_UNIQUE"_f1"
time $PRG $SAMPLE2_NOVEL_UNIQUE $M_MIN_NSR2 $M_MIN_NRS $M_MIN_NRD2 $SAMPLE2_NOVEL_UNIQUE"_f1"

file_size=(`du -cb $SAMPLE2_NOVEL_UNIQUE"_f1"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFiltermean3sd/./SXFiltermean3sd" 

echo ""
echo $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f1" $M_Mean3SD_Common $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f2"
time $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f1" $M_Mean3SD_Common $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f2"

file_size=(`du -cb $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f2"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_NOVEL_UNIQUE"_f1" $M_Mean3SD_1 $SAMPLE1_NOVEL_UNIQUE"_f2"
time $PRG $SAMPLE1_NOVEL_UNIQUE"_f1" $M_Mean3SD_1 $SAMPLE1_NOVEL_UNIQUE"_f2"

file_size=(`du -cb $SAMPLE1_NOVEL_UNIQUE"_f2"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_NOVEL_UNIQUE"_f1" $M_Mean3SD_2 $SAMPLE2_NOVEL_UNIQUE"_f2"
time $PRG $SAMPLE2_NOVEL_UNIQUE"_f1" $M_Mean3SD_2 $SAMPLE2_NOVEL_UNIQUE"_f2"

file_size=(`du -cb $SAMPLE2_NOVEL_UNIQUE"_f2"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi


#------------------------------------------------------------------------------------------
PRG=$M_APPDIR"/SXFilterRepeat/./SXFilterRepeat"

echo ""
echo $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f2" $M_MASKDENSDIR $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f3"
time $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f2" $M_MASKDENSDIR $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f3"

file_size=(`du -cb $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f3"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_NOVEL_UNIQUE"_f2" $M_MASKDENSDIR $SAMPLE1_NOVEL_UNIQUE"_f3"
time $PRG $SAMPLE1_NOVEL_UNIQUE"_f2" $M_MASKDENSDIR $SAMPLE1_NOVEL_UNIQUE"_f3"

file_size=(`du -cb $SAMPLE1_NOVEL_UNIQUE"_f3"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_NOVEL_UNIQUE"_f2" $M_MASKDENSDIR $SAMPLE2_NOVEL_UNIQUE"_f3"
time $PRG $SAMPLE2_NOVEL_UNIQUE"_f2" $M_MASKDENSDIR $SAMPLE2_NOVEL_UNIQUE"_f3"

file_size=(`du -cb $SAMPLE2_NOVEL_UNIQUE"_f3"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi


#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXFilterFinalSRGT/./SXFilterFinalSRGT"

SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED=$SAMPLE1_SAMPLE2_NOVEL_COMMON"_f"$M_MAX_SR
SAMPLE1_NOVEL_UNIQUE_FILTERED=$SAMPLE1_NOVEL_UNIQUE"_f"$M_MAX_SR
SAMPLE2_NOVEL_UNIQUE_FILTERED=$SAMPLE2_NOVEL_UNIQUE"_f"$M_MAX_SR

echo ""
echo $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f3" $SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED $M_MAX_SR
time $PRG $SAMPLE1_SAMPLE2_NOVEL_COMMON"_f3" $SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED $M_MAX_SR

file_size=(`du -cb $SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE1_NOVEL_UNIQUE"_f3" $SAMPLE1_NOVEL_UNIQUE_FILTERED $M_MAX_SR
time $PRG $SAMPLE1_NOVEL_UNIQUE"_f3" $SAMPLE1_NOVEL_UNIQUE_FILTERED $M_MAX_SR

file_size=(`du -cb $SAMPLE1_NOVEL_UNIQUE_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo $PRG $SAMPLE2_NOVEL_UNIQUE"_f3" $SAMPLE2_NOVEL_UNIQUE_FILTERED $M_MAX_SR
time $PRG $SAMPLE2_NOVEL_UNIQUE"_f3" $SAMPLE2_NOVEL_UNIQUE_FILTERED $M_MAX_SR

file_size=(`du -cb $SAMPLE2_NOVEL_UNIQUE_FILTERED`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FILTER_LIST_1=$M_SAMPLEID_1"_"$SNAME"_Sample"

eval cat $SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED $SAMPLE1_KNOWN_UNIQUE_FILTERED $SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED $SAMPLE1_NOVEL_UNIQUE_FILTERED ">" $FILTER_LIST_1

file_size=(`du -cb $FILTER_LIST_1`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

FILTER_LIST_2=$M_SAMPLEID_2"_"$SNAME"_Sample"

eval cat $SAMPLE1_SAMPLE2_KNOWN_COMMON_FILTERED $SAMPLE2_KNOWN_UNIQUE_FILTERED  $SAMPLE1_SAMPLE2_NOVEL_COMMON_FILTERED $SAMPLE2_NOVEL_UNIQUE_FILTERED ">" $FILTER_LIST_2

file_size=(`du -cb $FILTER_LIST_2`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

#if [ "$M_TYPE" = "S" ]; then
   #M_FILTER_LIST=$FILTER_LIST_1 #".lst"

#   PRG=$M_APPDIR"/SXCalculateRSGraph2d/./SXCalculateRSGraph2d"

#   echo ""    
#   echo $PRG $FILTER_LIST_1 $FILTER_LIST_1".Graph"
#   time $PRG $FILTER_LIST_1 $FILTER_LIST_1".Graph"

#   file_size=(`du -cb $FILTER_LIST_1".Graph"`)
#   if [ $file_size -eq 0 ]; then
#      exit 1
#   fi

#   echo ""    
#   echo $PRG $FILTER_LIST_2 $FILTER_LIST_2".Graph"
#   time $PRG $FILTER_LIST_2 $FILTER_LIST_2".Graph"

#   file_size=(`du -cb $FILTER_LIST_2".Graph"`)
#   if [ $file_size -eq 0 ]; then
#      exit 1
#   fi
#fi


