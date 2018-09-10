#!/bin/sh

SAMPLE_ID=$1
SAMPLE_NO=$2

FULL_LIST=$SAMPLE_ID"_snv_close_vicinity"

SNP_AQS_LIST=$SAMPLE_ID"_snp_Sample"

INS_MERGEPOLY_LIST=$SAMPLE_ID"_ins_c1_aqs_110118.lst"
DEL_MERGEPOLY_LIST=$SAMPLE_ID"_del_c1_aqs_110118.lst"

INS_LIST=$SAMPLE_ID"_snv_close_vicinity_INS"
DEL_LIST=$SAMPLE_ID"_snv_close_vicinity_DEL"

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtCloseVicinitySNV"

echo $PRG "../"$SNP_AQS_LIST $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST  
$PRG "../"$SNP_AQS_LIST $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST

#------------------------------------------------------------------------------------------------

echo "grep \"D-0\"" FULL_LIST ">" $INS_LIST
grep "D-0" $FULL_LIST > $INS_LIST

echo "grep \"I-0\"" FULL_LIST ">" $DEL_LIST
grep "I-0" $FULL_LIST > $DEL_LIST

#------------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtPolybaseRecs"

echo $PRG "../""PolyINS_"$SAMPLE_ID".lst" $INS_LIST "I" $INS_LIST".poly" 
time $PRG "../""PolyINS_"$SAMPLE_ID".lst" $INS_LIST "I" $INS_LIST".poly"

echo $PRG "../""PolyDEL_"$SAMPLE_ID".lst" $DEL_LIST "D" $DEL_LIST".poly" 
time $PRG "../""PolyDEL_"$SAMPLE_ID".lst" $DEL_LIST "D" $DEL_LIST".poly"

#------------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXFilterSNPPoly"

CURR_DATE=$(date '+%y%m%d')
#SNP_AQS_LIST=$SAMPLE_ID"_snp_c1_aqs_"CURR_DATE".lst"
#SNP_AQS_LIST=$SAMPLE_ID"_snp_c1_aqs_110117.lst"

echo $PRG $INS_LIST".poly" "../"$SNP_AQS_LIST $SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INS_filtered"
time $PRG $INS_LIST".poly" "../"$SNP_AQS_LIST $SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INS_filtered"

echo $PRG $DEL_LIST".poly" $SNP_AQS_LIST"_INS" $SAMPLE_ID"_snp_c1_aqs_"$CURR_DATE".lst" $SNP_AQS_LIST"_DEL_filtered"
time $PRG $DEL_LIST".poly" $SNP_AQS_LIST"_INS" $SAMPLE_ID"_snp_c1_aqs_"$CURR_DATE".lst" $SNP_AQS_LIST"_DEL_filtered"

#------------------------------------------------------------------------------------------------

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXMergeSNPSample2INDEL"

#INS_LIST=$SAMPLE_ID"_ins_c1_aqs_"CURR_DATE".lst"

echo $PRG $SNP_AQS_LIST"_INS_filtered" $SAMPLE_NO $INS_LIST".poly" $INS_MERGEPOLY_LIST "I" $SAMPLE_ID"_ins_c1_aqs_"$CURR_DATE".lst"
time $PRG $SNP_AQS_LIST"_INS_filtered" $SAMPLE_NO $INS_LIST".poly" $INS_MERGEPOLY_LIST "I" $SAMPLE_ID"_ins_c1_aqs_"$CURR_DATE".lst"

echo $PRG $SNP_AQS_LIST"_DEL_filtered" $SAMPLE_NO $DEL_LIST".poly" $DEL_MERGEPOLY_LIST "D" $SAMPLE_ID"_del_c1_aqs_"$CURR_DATE".lst"
time $PRG $SNP_AQS_LIST"_DEL_filtered" $SAMPLE_NO $DEL_LIST".poly" $DEL_MERGEPOLY_LIST "D" $SAMPLE_ID"_del_c1_aqs_"$CURR_DATE".lst"






