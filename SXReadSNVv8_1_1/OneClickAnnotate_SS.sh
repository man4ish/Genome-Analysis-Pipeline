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
#if [ "$M_TYPE" = "" ]; then TYPE=$5; else TYPE=$M_TYPE; fi
#if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$6; else SAMPLEID=$M_SAMPLEID; fi
#if [ "$M_GENOME" = "" ]; then GENOME_VER=$7; else GENOME_VER=$M_GENOME; fi

MIN_PSR=1
MAX_PSR=500
TYPE=$3 
SAMPLEID=$4
GENOME_VER=$M_GENOME

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

FBINNED=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
FEXCEPT=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
FANNOTATE=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVv6/./SXReadSNV"

echo $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $TYPE $SAMPLEID $GENOME_VER
time $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $TYPE $SAMPLEID $GENOME_VER

if [ "$TYPE" = "S" ]; then
   if [ "$GENOME_VER" = "37.1" ]; then  
      DB=$DB_37_1 
   else
      DB=$DB_36_3
   fi

   echo $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"
   time $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"
fi

echo "Done...!!!"
