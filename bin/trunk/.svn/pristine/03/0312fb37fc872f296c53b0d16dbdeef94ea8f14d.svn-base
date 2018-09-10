#!/bin/sh

APPDIR="/archive/project/oneclick/syhwah"

#if [ "$1" != "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
#if [ "$2" != "" ]; then MIN_PSR=$2; else MIN_PSR=$M_MIN_PSR; fi
#if [ "$3" != "" ]; then MAX_PSR=$3; else MAX_PSR=$M_MAX_PSR; fi
#if [ "$4" != "" ]; then TYPE=$4; else TYPE=$M_TYPE; fi
#if [ "$5" != "" ]; then SAMPLEID=$5; else M_SAMPLEID=$M_SAMPLEID; fi
#if [ "$6" != "" ]; then GENOME_VER=$6; else GENOME_VER=$M_GENOME; fi


if [ "$M_FILTER_LIST" = "" ]; then FILE_INPUT=$1; else FILE_INPUT=$M_FILTER_LIST; fi
if [ "$M_MIN_PSR" = "" ]; then MIN_PSR=$2; else MIN_PSR=$M_MIN_PSR; fi
if [ "$M_MAX_PSR" = "" ]; then MAX_PSR=$3; else MAX_PSR=$M_MAX_PSR; fi
if [ "$M_TYPE" = "" ]; then TYPE=$4; else TYPE=$M_TYPE; fi
if [ "$M_SAMPLEID" = "" ]; then SAMPLEID=$5; else M_SAMPLEID=$M_SAMPLEID; fi
if [ "$M_GENOME" = "" ]; then GENOME_VER=$6; else GENOME_VER=$M_GENOME; fi

REFDIR=$APPDIR"/Genome_"$GENOME_VER
FGENE=$REFDIR"/Genedb_"$GENOME_VER
FEXON=$REFDIR"/exon_"$GENOME_VER"_ccds" #CCDS_Exon
FCNV=$REFDIR"/cnv_Variation_"$GENOME_VER"_ncbi"

PRG=$APPDIR"/SXReadSNVv4/./SXReadSNV"

echo $PRG "-f" $FILE_INPUT "-gL" $FGENE $FEXON $FCNV $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID
time $PRG -f $FILE_INPUT -gL $FGENE $FEXON $FCNV $MIN_PSR $MAX_PSR $M_TYPE $M_SAMPLEID

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
