#!/bin/sh

FULL_LIST=$M_SAMPLEID"_snv_close_vicinity"
SNP_SAMPLE=$M_SAMPLEID"_snp_Sample"

CURR_DATE=$(date '+%y%m%d')

INS_MERGEPOLY_LIST=$M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst.merge"

#file_size=(`du -cb $INS_MERGEPOLY_LIST`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

DEL_MERGEPOLY_LIST=$M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst.merge"

#file_size=(`du -cb $DEL_MERGEPOLY_LIST`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

INS_LIST=$M_SAMPLEID"_snv_close_vicinity_INS"
DEL_LIST=$M_SAMPLEID"_snv_close_vicinity_DEL"

#-----------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/CloseVicinity/./SXExtCloseVicinitySNV"

echo $PRG $SNP_SAMPLE $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST  
$PRG $SNP_SAMPLE $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST

#file_size=(`du -cb $FULL_LIST`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

echo "grep \"D-0\"" FULL_LIST ">" $INS_LIST
grep "D-0" $FULL_LIST > $INS_LIST

echo "grep \"I-0\"" FULL_LIST ">" $DEL_LIST
grep "I-0" $FULL_LIST > $DEL_LIST

#-----------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtPolybaseRecs"

echo $PRG "PolyINS_"$M_SAMPLEID".lst" $INS_LIST "I" $INS_LIST".poly" 
time $PRG "PolyINS_"$M_SAMPLEID".lst" $INS_LIST "I" $INS_LIST".poly"

#file_size=(`du -cb $INS_LIST".poly"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

echo $PRG "PolyDEL_"$M_SAMPLEID".lst" $DEL_LIST "D" $DEL_LIST".poly" 
time $PRG "PolyDEL_"$M_SAMPLEID".lst" $DEL_LIST "D" $DEL_LIST".poly"

#file_size=(`du -cb $DEL_LIST".poly"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#-----------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXFilterSNPPoly"

echo $PRG $INS_LIST".poly" $SNP_SAMPLE $SNP_SAMPLE"_INS" $SNP_SAMPLE"_INS_filtered"
time $PRG $INS_LIST".poly" $SNP_SAMPLE $SNP_SAMPLE"_INS" $SNP_SAMPLE"_INS_filtered"

#file_size=(`du -cb $SNP_SAMPLE"_INS_filtered"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

mv $SNP_SAMPLE $SNP_SAMPLE".org"

echo $PRG $DEL_LIST".poly" $SNP_SAMPLE"_INS" $SNP_SAMPLE $SNP_SAMPLE"_DEL_filtered"
time $PRG $DEL_LIST".poly" $SNP_SAMPLE"_INS" $SNP_SAMPLE $SNP_SAMPLE"_DEL_filtered"

#file_size=(`du -cb $SNP_SAMPLE"_DEL_filtered"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#-----------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXMergeSNPSample2INDEL"

echo $PRG $SNP_SAMPLE"_INS_filtered" $M_SAMPLENO $INS_LIST".poly" $INS_MERGEPOLY_LIST "I" $M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst"
time $PRG $SNP_SAMPLE"_INS_filtered" $M_SAMPLENO $INS_LIST".poly" $INS_MERGEPOLY_LIST "I" $M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst"

#file_size=(`du -cb $M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

echo $PRG $SNP_SAMPLE"_DEL_filtered" $M_SAMPLENO $DEL_LIST".poly" $DEL_MERGEPOLY_LIST "D" $M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst"
time $PRG $SNP_SAMPLE"_DEL_filtered" $M_SAMPLENO $DEL_LIST".poly" $DEL_MERGEPOLY_LIST "D" $M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst"

#file_size=(`du -cb $M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi




