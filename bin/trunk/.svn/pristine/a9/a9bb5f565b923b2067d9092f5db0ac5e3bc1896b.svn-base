#!/bin/sh

case "$M_TYPE" in
   "I")
       SNAME="ins"
       ;;
   "D")
       SNAME="del"
       ;;
esac

FILE_INPUT=$M_SAMPLEID"_"$SNAME"_Sample"
READSDIR=$M_APPDIR"/GS000117/37_1_"$M_SAMPLEID"/raw_read"
PERFECTSDIR=$READSDIR
REPTDENDIR=$READSDIR

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendPerfectDenPoly/./SXAppendPDenMinus2Poly_SingleSample"

echo ""
echo $PRG $FILE_INPUT "PDens_"$SNAME"_"$M_SAMPLEID".lst" $M_SAMPLENO $FILE_INPUT"_perfect"
time $PRG $FILE_INPUT "PDens_"$SNAME"_"$M_SAMPLEID".lst" $M_SAMPLENO $FILE_INPUT"_perfect"

file_size=(`du -cb $FILE_INPUT"_perfect"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/Java_SXAppendZygo/CGZygosity.java"

#PRG="java CGZygosity"

#echo ""
#echo $PRG $FILE_INPUT"_perfect" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"
#time $PRG $FILE_INPUT"_perfect" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"


case "$M_TYPE" in
   "I")
       PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityINS"
       ;;
   "D")
       PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityDEL"
       ;;
esac

echo ""
echo $PRG $FILE_INPUT"_perfect" $M_SAMPLENO $FILE_INPUT"_Zygo"
time $PRG $FILE_INPUT"_perfect" $M_SAMPLENO $FILE_INPUT"_Zygo"


file_size=(`du -cb $FILE_INPUT"_Zygo"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintCNV"

echo ""
echo $PRG $FILE_INPUT"_Zygo" $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $FILE_INPUT"_ZygoCNV"
time $PRG $FILE_INPUT"_Zygo" $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $FILE_INPUT"_ZygoCNV"

file_size=(`du -cb $FILE_INPUT"_ZygoCNV"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FBINNED=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVvPoly2_dbSNP131/./SXReadSNV"

echo ""
echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME

#------------------------------------------------------------------------------------------

CURR_DATE=$(date '+%y%m%d')

if [ $M_TYPE="I" ]; then
   PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreINS"
else
   PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreDEL"
fi

IN_FILE=$M_SAMPLEID"_"$SNAME"_c1_"$CURR_DATE".lst"
OUT_FILE=$IN_FILE".score"

echo ""
echo $PRG $IN_FILE $M_MEAN $OUT_FILE
time $PRG $IN_FILE $M_MEAN $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

case "$M_TYPE" in
   "I")
       PRG=$M_APPDIR"/SXPrintStats/SXPrintStatsINS"
       ;;
   "D")
       PRG=$M_APPDIR"/SXPrintStats/SXPrintStatsDEL"
       ;;
esac
   
echo $PRG $OUT_FILE "1" $OUT_FILE".stats.rpt"
time $PRG $OUT_FILE "1" $OUT_FILE".stats.rpt"
   
file_size=(`du -cb $OUT_FILE".stats.rpt"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendDrugInfo/./SXAppendDrugInfo"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".druginfo"

echo ""
echo $PRG "dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE
time $PRG $M_APPDIR"/Genome_37.1/dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo "Done...!!!"
