#!/bin/sh

#M_GENOME="36.3"
M_APPDIR="/archive/project/oneclick/syhwah"

#if [ "$M_FILTER_LIST" = "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
#if [ "$M_ZYGO_CUTOFF" = "" ]; then ZYGO_CUTOFF=$2; else ZYGO_CUTOFF=$M_ZYGO_CUTOFF; fi

#PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV_SingleSample" 

#echo $PRG $FILE_INPUT $ZYGO_CUTOFF $FILE_INPUT".ZygoCNV"
#time $PRG $FILE_INPUT $ZYGO_CUTOFF $FILE_INPUT".ZygoCNV"

#if [ "$M_TYPE" = "" ]; then TYPE=$3; else TYPE=$M_TYPE; fi
#if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$4; else SAMPLEID=$M_SAMPLEID; fi
#if [ "$M_GENOME" = "" ]; then GENOME_VER=$7; else GENOME_VER=$M_GENOME; fi

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

PRG="java -jar "$M_APPDIR"/Java_SXAppendZygo/CGZygosity_SingleSample.jar"

echo ""
echo $PRG $FILE_INPUT"_perfect" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"
time $PRG $FILE_INPUT"_perfect" $M_ZYGO_MIN_SR $M_SAMPLENO $M_DIPLOID_PROB $FILE_INPUT"_Zygo"

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


MIN_PSR=1
MAX_PSR=500

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

if [ $M_GENOME = "36.3" ]; then
   FBINNED=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
   FEXCEPT=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
   FANNOTATE=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt
else
   FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
   FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
   FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt
fi

PRG=$M_APPDIR"/SXReadSNVvSIMPolyPDens2/./SXReadSNVSim"

echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID $M_GENOME
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID $M_GENOME


#case "$M_TYPE" in
#   "I")
#       PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreINS_SingleSampleSIM"  
#       ;;
#   "D")
#       PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreDEL_SingleSampleSIM"
#       ;;
#esac

PRG=$M_APPDIR"/SXAppendScore/./SXAppendScore_SingleSampleSIMCG"

CURR_DATE=$(date '+%y%m%d')
IN_FILE=$M_SAMPLEID"_"$SNAME"_c1_"$CURR_DATE".lst"
OUT_FILE=$IN_FILE".score"

echo $PRG $IN_FILE $M_MEAN $OUT_FILE "-f"
time $PRG $IN_FILE $M_MEAN $OUT_FILE -f

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------


echo "Done...!!!"
