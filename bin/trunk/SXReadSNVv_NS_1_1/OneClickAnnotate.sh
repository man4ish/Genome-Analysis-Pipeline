#!/bin/sh

M_GENOME="36.3"
M_CG_SAMPLEID="NA19240"
M_APPDIR="/archive/project/oneclick/syhwah"

if [ $M_GENOME = "37.1" ]; then
   READSDIR=$M_APPDIR"/GS000117/37_1_"$M_CG_SAMPLEID"/raw_read"
elif [ $M_GENOME = "36.3" ]; then
   READSDIR=$M_APPDIR"/GS000117/"$M_CG_SAMPLEID"/raw_read"
else
   READSDIR=$M_APPDIR"/Illumina/Harvard/run_2/raw_read"
fi

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

#if [ "$M_FILTER_LIST" = "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
#if [ "$M_ZYGO_CUTOFF" = "" ]; then $M_ZYGO_CUTOFF=$2; fi

FILE_INPUT=$1
ZYGO_CUTOFF=$2
MEANS=$3
M_TYPE=$4
M_SAMPLEID=$5

HET_MSR=4
HOM_MSR=5

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV" 

echo $FILE_INPUT $ZYGO_CUTOFF $HET_MSR $HOM_MSR $MEANS $READSDIR $FILE_INPUT".ZygoCNV"
time $FILE_INPUT $ZYGO_CUTOFF $HET_MSR $HOM_MSR $MEANS $READSDIR $FILE_INPUT".ZygoCNV"

#if [ "$M_MIN_PSR" = "" ]; then MIN_PSR=$3; else MIN_PSR=$M_MIN_PSR; fi
#if [ "$M_MAX_PSR" = "" ]; then MAX_PSR=$4; else MAX_PSR=$M_MAX_PSR; fi
#if [ "$M_TYPE" = "" ]; then TYPE=$5; else TYPE=$M_TYPE; fi
#if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$6; else M_SAMPLEID=$M_SAMPLEID; fi
#if [ "$M_GENOME" = "" ]; then GENOME_VER=$7; else GENOME_VER=$M_GENOME; fi

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$APPDIR"/SXReadSNVv6/./SXReadSNV"

#echo $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID
#time $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID

echo $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID
time $PRG "-f" $FILE_INPUT".ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID

if [ "$TYPE" = "S" ]; then
   if [ "$GENOME_VER" = "37.1" ]; then  
      DB=$DB_37_1 
   else
      DB=$DB_36_3
   fi

   echo $APPDIR"/SXSynoMap/./SXSynoMap" $DB $M_SAMPLEID".syn" $M_SAMPLEID".fna"
   time $APPDIR"/SXSynoMap/./SXSynoMap" $DB $M_SAMPLEID".syn" $M_SAMPLEID".fna"
fi

echo "Done...!!!"
