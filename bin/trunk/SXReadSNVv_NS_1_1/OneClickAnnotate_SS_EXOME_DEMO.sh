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
REPTDENDIR=$M_READSDIR

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendPerfectDen/./SXAppendReadDenMinus2"

echo ""
echo $PRG $FILE_INPUT $M_READSDIR $M_TYPE $M_SAMPLENO $FILE_INPUT"_ReadDen_Minus2"
time $PRG $FILE_INPUT $M_READSDIR $M_TYPE $M_SAMPLENO $FILE_INPUT"_ReadDen_Minus2"

file_size=(`du -cb $FILE_INPUT"_ReadDen_Minus2"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/Java_SXAppendZygo/CGZygosity.java"

PRG="java CGZygosity_SingleSample"

echo ""
echo $PRG $FILE_INPUT"_ReadDen_Minus2" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"
time $PRG $FILE_INPUT"_ReadDen_Minus2" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"

file_size=(`du -cb $FILE_INPUT"_Zygo"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintCNVNS"

echo ""
echo $PRG $FILE_INPUT"_Zygo 0" $M_READSDIR $REPTDENDIR $FILE_INPUT"_ZygoCNV"
time $PRG $FILE_INPUT"_Zygo" 0 $M_READSDIR $REPTDENDIR $FILE_INPUT"_ZygoCNV"

file_size=(`du -cb $FILE_INPUT"_ZygoCNV"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVv_NS_1_1/./SXReadSNV"

echo ""
echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME

#------------------------------------------------------------------------------------------

CURR_DATE=$(date '+%y%m%d')
PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreNS"
IN_FILE=$M_SAMPLEID"_"$SNAME"_c1_"$CURR_DATE".lst"
OUT_FILE=$IN_FILE".score"

echo ""
echo $PRG $IN_FILE $OUT_FILE
time $PRG $IN_FILE $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendDrugInfo/./SXAppendDrugInfo"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".druginfo"

echo ""
echo $PRG  $M_APPDIR"dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE
time $PRG $M_APPDIR"/Genome_37.1/dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

echo ""
echo "Done...!!!"
