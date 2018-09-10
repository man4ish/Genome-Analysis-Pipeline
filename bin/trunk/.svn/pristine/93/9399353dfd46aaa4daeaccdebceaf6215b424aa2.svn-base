#!/bin/sh

FULL_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_snv_close_vicinity"
SNP_SAMPLE=$M_OUTPUT_DIR"/"$M_SAMPLEID"_snp_Sample"

CURR_DATE=$(date '+%y%m%d')

INS_MERGEPOLY_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst.merge"

file_size=(`du -cb $INS_MERGEPOLY_LIST`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

DEL_MERGEPOLY_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst.merge"

file_size=(`du -cb $DEL_MERGEPOLY_LIST`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

INS_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_snv_close_vicinity_INS"
DEL_LIST=$M_OUTPUT_DIR"/"$M_SAMPLEID"_snv_close_vicinity_DEL"

#-----------------------------------------------------------------------------------------------

PRG="SXExtCloseVicinitySNV"

eval /usr/bin/time --verbose $PRG $SNP_SAMPLE $INS_MERGEPOLY_LIST $DEL_MERGEPOLY_LIST $FULL_LIST \
                             >> $M_LOG_FILE 2>> $M_LOG_FILE

#file_size=(`du -cb $FULL_LIST`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

eval /usr/bin/time --verbose grep "D-0" $FULL_LIST ">" $INS_LIST 2>> $M_LOG_FILE
eval /usr/bin/time --verbose grep "I-0" $FULL_LIST ">" $DEL_LIST 2>> $M_LOG_FILE

#-----------------------------------------------------------------------------------------------

PRG="SXExtPolybaseRecs"

eval /usr/bin/time --verbose $PRG $M_OUTPUT_DIR"/PolyINS_"$M_SAMPLEID".lst" $INS_LIST "I" $INS_LIST".poly" >> $M_LOG_FILE \
                             2>> $M_LOG_FILE

#file_size=(`du -cb $INS_LIST".poly"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

eval /usr/bin/time --verbose $PRG $M_OUTPUT_DIR"/PolyDEL_"$M_SAMPLEID".lst" $DEL_LIST "D" $DEL_LIST".poly" >> $M_LOG_FILE \
                             2>> $M_LOG_FILE

#file_size=(`du -cb $DEL_LIST".poly"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#-----------------------------------------------------------------------------------------------

PRG="SXFilterSNPPoly"

eval /usr/bin/time --verbose $PRG $INS_LIST".poly" $SNP_SAMPLE $SNP_SAMPLE"_INS" $SNP_SAMPLE"_INS_filtered" \
                             >> $M_LOG_FILE 2>> $M_LOG_FILE

#file_size=(`du -cb $SNP_SAMPLE"_INS_filtered"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

mv $SNP_SAMPLE $SNP_SAMPLE".org"

eval /usr/bin/time --verbose $PRG $DEL_LIST".poly" $SNP_SAMPLE"_INS" $SNP_SAMPLE $SNP_SAMPLE"_DEL_filtered" \
                             >> $M_LOG_FILE 2>> $M_LOG_FILE

#file_size=(`du -cb $SNP_SAMPLE"_DEL_filtered"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#-----------------------------------------------------------------------------------------------

PRG="SXMergeSNPSample2INDEL_SingleSample"

eval /usr/bin/time --verbose $PRG $SNP_SAMPLE"_INS_filtered" $M_SAMPLENO $INS_LIST".poly" $INS_MERGEPOLY_LIST "I" \
                             $M_OUTPUT_DIR"/"$M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst" >> $M_LOG_FILE 2>> $M_LOG_FILE

#file_size=(`du -cb $M_OUTPUT_DIR"/"$M_SAMPLEID"_ins_c1_aqs_"$CURR_DATE".lst"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

eval /usr/bin/time --verbose $PRG $SNP_SAMPLE"_DEL_filtered" $M_SAMPLENO $DEL_LIST".poly" $DEL_MERGEPOLY_LIST "D" \
                            $M_OUTPUT_DIR"/"$M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst" >> $M_LOG_FILE 2>> $M_LOG_FILE

#file_size=(`du -cb $M_OUTPUT_DIR"/"$M_SAMPLEID"_del_c1_aqs_"$CURR_DATE".lst"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi




