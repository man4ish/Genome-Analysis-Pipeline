#!/bin/sh

M_GENOME="36.3"
M_APPDIR="/archive/project/oneclick/syhwah"

if [ "$M_FILTER_LIST" = "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
if [ "$M_ZYGO_CUTOFF" = "" ]; then ZYGO_CUTOFF=$2; else ZYGO_CUTOFF=$M_ZYGO_CUTOFF; fi

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV_SingleSample" 

echo $PRG $FILE_INPUT $ZYGO_CUTOFF $FILE_INPUT".ZygoCNV"
time $PRG $FILE_INPUT $ZYGO_CUTOFF $FILE_INPUT".ZygoCNV"

#if [ "$M_MIN_PSR" = "" ]; then MIN_PSR=$3; else MIN_PSR=$M_MIN_PSR; fi
#if [ "$M_MAX_PSR" = "" ]; then MAX_PSR=$4; else MAX_PSR=$M_MAX_PSR; fi
if [ "$M_TYPE" = "" ]; then TYPE=$3; else TYPE=$M_TYPE; fi
if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$4; else SAMPLEID=$M_SAMPLEID; fi
if [ "$M_GENOME" = "" ]; then GENOME_VER=$7; else GENOME_VER=$M_GENOME; fi

MIN_PSR=1
MAX_PSR=500
#TYPE=$3 
#SAMPLEID=$4
#GENOME_VER=$M_GENOME

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

FBINNED=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
FEXCEPT=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
FANNOTATE=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVv8_1_2/./SXReadSNV"

echo $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $TYPE $SAMPLEID $GENOME_VER
time $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $TYPE $SAMPLEID $GENOME_VER

CURR_DATE=$(date '+%y%m%d')
SNAME="snp"

IN_FILE=$M_SAMPLEID"_"$SNAME"_c1_"$CURR_DATE

if [ "$GENOME_VER" = "37.1" ]; then  
     DB=$DB_37_1 
else
     DB=$DB_36_3
fi

echo $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"
time $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"

#------------------------------------------------------------------------------------------
#synasearch 
#
#echo $SYNASEARCH -np -f $M_SAMPLEID".fna" -db $DB > $M_SAMPLEID".sso"
#time eval $SYNASEARCH -np -f $M_SAMPLEID".fna" -db $DB > $M_SAMPLEID".sso"
#
#file_size=(`du -cb $M_SAMPLEID".sso"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/SSNpParser/./SSNpParser"

#OUT_FILE=$IN_FILE #".lst"

#echo $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID  
#time $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID

#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#     exit 1
#fi

#IN_FILE=$OUT_FILE
#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreSNP_SS_wo_PVInfo"
#
#OUT_FILE=$IN_FILE".score"
#
#echo $PRG $IN_FILE $OUT_FILE
#time $PRG $IN_FILE $OUT_FILE
#
#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------


echo "Done...!!!"
