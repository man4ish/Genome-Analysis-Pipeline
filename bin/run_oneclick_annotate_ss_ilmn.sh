#!/bin/sh

case "$M_TYPE" in
   "I")
       SNAME="ins"
       ;;
   "D")
       SNAME="del"
       ;;
esac

REPTDENDIR=$M_READSDIR
OUT_FILE=$M_FILTER_LIST"_perfect"

#------------------------------------------------------------------------------------------

PRG="SXAppendPDenMinus2Poly_SingleSample"

OUT_FILE=$M_FILTER_LIST"_perfect"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $M_FILTER_LIST $M_OUTPUT_DIR"/PDens_"$SNAME"_"$M_SAMPLEID".lst" \
                             $M_SAMPLENO $OUT_FILE >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG="SXPZygoINDEL_SingleSample"

IN_FILE=$OUT_FILE
OUT_FILE=$M_FILTER_LIST"_Zygo"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $IN_FILE $M_MEAN $OUT_FILE >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG="SXPrintCNVNS"

IN_FILE=$OUT_FILE
OUT_FILE=$M_FILTER_LIST"_ZygoCNV"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $IN_FILE $M_SAMPLEID 0 $M_READSDIR $REPTDENDIR $OUT_FILE \
                             >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

CURR_DATE=$(date '+%y%m%d')

case "$M_TYPE" in
   "I")        
       #PRG="SXAppendScoreINS_IlmnSim_B4Annotation"
       PRG="SXAppendScoreINS_Ilmn_B4Annotation"
       ;;
   "D")
       #PRG="SXAppendScoreDEL_IlmnSim_B4Annotation"
       PRG="SXAppendScoreDEL_Ilmn_B4Annotation" 
       ;;
esac

IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".score"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $IN_FILE $M_MEAN $OUT_FILE >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

INFILE=$OUT_FILE
FILTER_LIST=$INFILE".filtered"

PRG="SXFilterByCScore"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $INFILE $M_MIN_CSCORE $FILTER_LIST  >> $M_LOG_FILE 2>> $M_LOG_FILE

#------------------------------------------------------------------------------------------

if [ $M_GENOME = "37.1" ]; then
   FBINNED=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
   FEXCEPT=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
   FANNOTATE=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt
   DB=$DB_37_1
else
   FBINNED=$M_REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
   FEXCEPT=$M_REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
   FANNOTATE=$M_REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt
   DB=$DB_36_3
fi

PRG="SXReadSNV_Ilm"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG "-f" $FILTER_LIST "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE \
                             $M_SAMPLEID $M_GENOME $M_FINALOUT_DIR $M_PROJ_SAMPLEID >> $M_LOG_FILE 2>> $M_LOG_FILE

#------------------------------------------------------------------------------------------

CURR_DATE=$(date '+%y%m%d')

case "$M_TYPE" in
   "I")
       #PRG="SXAppendScoreINS_IlmnSim"
       PRG="SXAppendScoreINS_Ilmn"
       ;;
   "D")
       #PRG="SXAppendScoreDEL_IlmnSim"
       PRG="SXAppendScoreDEL_Ilmn" 
       ;;
esac

IN_FILE=$M_FINALOUT_DIR"/"$M_SAMPLEID"_"$SNAME"_c1_"$CURR_DATE".lst"
OUT_FILE=$IN_FILE".score"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG  $IN_FILE $M_MEAN $OUT_FILE >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

case "$M_TYPE" in
   "I")
       PRG="SXPrintStatsINS"
       ;;
   "D")
       PRG="SXPrintStatsDEL"
       ;;
esac

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $OUT_FILE "1" $OUT_FILE".stats.rpt" >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE".stats.rpt"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
#Update Summary Report

PRG="Avg_Confi_Score"

SUMMARY_RPT=$M_FINALOUT_DIR"/"$M_SAMPLEID"_"$SNAME"_summary_c1_"$CURR_DATE".rpt"
SUMMARY_ORI=$M_OUTPUT_DIR"/"$M_SAMPLEID"_"$SNAME"_summary_c1_"$CURR_DATE".rpt.ori"

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose mv $SUMMARY_RPT $SUMMARY_ORI 2>> $M_LOG_FILE

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $OUT_FILE $SUMMARY_ORI $SUMMARY_RPT >> $M_LOG_FILE 2>> $M_LOG_FILE                             

#------------------------------------------------------------------------------------------

PRG="SXAppendDrugInfo"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".druginfo"

if [ $M_GENOME = "37.1" ]; then
   DRUGINFO_FILE="dbSNP131_37.1_hg19_17August_Synamatix.txt" 
else
   DRUGINFO_FILE="dbSNP129_36.3_20Sept_Synamatix.txt" 
fi

echo "" >> $M_LOG_FILE
eval /usr/bin/time --verbose $PRG $M_REFDIR"/"$DRUGINFO_FILE $IN_FILE $OUT_FILE \
                             >> $M_LOG_FILE 2>> $M_LOG_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#echo ""
#echo "Done...!!!" >> $M_LOG_FILE
