#!/bin/sh

M_GENOME="37.1"
M_APPDIR="/archive/project/oneclick/syhwah"
REPTDENDIR=$M_READSDIR


if [ "$M_FILTER_LIST" = "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
if [ "$M_ZYGO_CUTOFF" = "" ]; then ZYGO_CUTOFF=$2; else ZYGO_CUTOFF=$M_ZYGO_CUTOFF; fi

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV_SingleSampleNS" 
ZERO=0
echo $PRG $FILE_INPUT $ZYGO_CUTOFF $ZERO $M_READSDIR $REPTDENDIR $FILE_INPUT".ZygoCNV"
time $PRG $FILE_INPUT $ZYGO_CUTOFF $ZERO $ZERO $ZERO $M_READSDIR $REPTDENDIR $FILE_INPUT".ZygoCNV"

if [ "$M_TYPE" = "" ]; then TYPE=$3; else TYPE=$M_TYPE; fi
if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$4; else SAMPLEID=$M_SAMPLEID; fi
if [ "$M_GENOME" = "" ]; then GENOME_VER=$7; else GENOME_VER=$M_GENOME; fi

MIN_PSR=1
MAX_PSR=500

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVv_NS_1_1/./SXReadSNV"

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

#echo $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"
#time $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"

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
