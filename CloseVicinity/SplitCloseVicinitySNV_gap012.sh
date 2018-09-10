#!/bin/sh

SAMPLE_ID=$1

FULL_LIST=$SAMPLE_ID"_snv_close_vicinity"

SNP_AQS_LIST=$SAMPLE_ID"_snp_c1_aqs_110117.lst"
INS_MERGEPOLY_LIST=$SAMPLE_ID"_ins_c1_aqs_110118.lst"
DEL_MERGEPOLY_LIST=$SAMPLE_ID"_del_c1_aqs_110118.lst"

INS_LIST=$SAMPLE_ID"_snv_close_vicinity_INS"
DEL_LIST=$SAMPLE_ID"_snv_close_vicinity_DEL"

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtCloseVicinitySNV"

echo $PRG "../"$SNP_AQS_LIST $INS_AQS_LIST $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST  
$PRG "../"$SNP_AQS_LIST $INS_AQS_LIST $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST

echo "grep \"D-0\"" FULL_LIST ">" $INS_LIST
#grep "D-0" $FULL_LIST > $INS_LIST

echo "grep \"I-0\"" FULL_LIST ">" $DEL_LIST
#grep "I-0" $FULL_LIST > $DEL_LIST

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtDataByGapSize" 

echo  $PRG $INS_LIST 0 I $INS_LIST"_gap0"
#time  $PRG $INS_LIST 0 "I" $INS_LIST"_gap0"

echo  $PRG $INS_LIST 1 I $INS_LIST"_gap1"
#time  $PRG $INS_LIST 1 "I" $INS_LIST"_gap1"

echo  $PRG $INS_LIST 2 I $INS_LIST"_gap2"
#time  $PRG $INS_LIST 2 "I" $INS_LIST"_gap2"

echo  $PRG $DEL_LIST 0 D $DEL_LIST"_gap0"
#time  $PRG $DEL_LIST 0 "D" $DEL_LIST"_gap0"

echo  $PRG $DEL_LIST 1 D $DEL_LIST"_gap1"
#time  $PRG $DEL_LIST 1 "D" $DEL_LIST"_gap1"

echo  $PRG $DEL_LIST 2 D $DEL_LIST"_gap2"
#time  $PRG $DEL_LIST 2 "D" $DEL_LIST"_gap2"

INS_LIST_GAPS=$INS_LIST"_gap01"

echo "cat " $INS_LIST"_gap0" $INS_LIST"_gap1 >" $INS_LIST_GAPS
#cat $INS_LIST"_gap0" $INS_LIST"_gap1 >" $INS_LIST_GAPS

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXExtPolybaseRecs"

echo $PRG "PolyINS_"$SAMPLE_ID".lst" $INS_LIST_GAPS "I" $INS_LIST_GAPS".poly" 
#time $PRG "PolyINS_"$SAMPLE_ID".lst" $INS_LIST_GAPS "I" $INS_LIST_GAPS".poly"

DEL_LIST_GAPS=$DEL_LIST"_gap01"

echo "cat " $DEL_LIST"_gap0" $DEL_LIST"_gap1 >" $DEL_LIST_GAPS
#cat $DEL_LIST"_gap0" $DEL_LIST"_gap1 >" $DEL_LIST_GAPS

echo $PRG "PolyDEL_"$SAMPLE_ID".lst" $DEL_LIST"_gap01 D" $DEL_LIST_GAPS".poly" 
#time $PRG "PolyDEL_"$SAMPLE_ID".lst" $DEL_LIST"_gap01 D" $DEL_LIST_GAPS".poly"

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXFilterSNPPoly"

CURR_DATE=$(date '+%y%m%d')
#SNP_AQS_LIST=$SAMPLE_ID"_snp_c1_aqs_"CURR_DATE".lst"
#SNP_AQS_LIST=$SAMPLE_ID"_snp_c1_aqs_110117.lst"

echo $PRG $INS_LIST_GAPS".poly" "../"$SNP_AQS_LIST $SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INS_filtered"
time $PRG $INS_LIST_GAPS".poly" "../"$SNP_AQS_LIST $SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INS_filtered"

echo $PRG $DEL_LIST_GAPS".poly" "../"$SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INSDEL" $SNP_AQS_LIST"_DEL_filtered"
time $PRG $DEL_LIST_GAPS".poly" "../"$SNP_AQS_LIST"_INS" $SNP_AQS_LIST"_INSDEL" $SNP_AQS_LIST"_DEL_filtered"

PRG="/archive/project/oneclick/syhwah/CloseVicinity/./SXMergeSNP2INDEL_Known"

#INS_LIST=$SAMPLE_ID"_ins_c1_aqs_"CURR_DATE".lst"
INS_LIST=$SAMPLE_ID"_ins_c1_aqs_110118.lst"

echo $PRG $SNP_LIST"_INS_filtered" $INS_LIST_GAPS".poly" $INS_LIST "I" $INS_LIST"_SNP"
time $PRG $SNP_LIST"_INS_filtered" $INS_LIST_GAPS".poly" $INS_LIST "I" $INS_LIST"_SNP"

#DEL_LIST=$SAMPLE_ID"_del_c1_aqs_"CURR_DATE".lst"
DEL_LIST=$SAMPLE_ID"_del_c1_aqs_110118.lst"

echo $PRG $SNP_LIST"_DEL_filtered" $DEL_LIST_GAPS".poly" $DEL_LIST "D" $DEL_LIST"_SNP"
time $PRG $SNP_LIST"_DEL_filtered" $DEL_LIST_GAPS".poly" $DEL_LIST "D" $DEL_LIST"_SNP"





