#!/bin/sh
export MALLOC_CHECK_=0
GENOME_VER="37.1"
REFDIR="/archive/project/oneclick/syhwah/Genome_"$GENOME_VER

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG="/archive/project/oneclick/syhwah/SXReadSNVv8_1_1/./SXReadSNV"
#PRG="./SXReadSNV"

echo $PRG "-f" $1 "-gL" $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER
time $PRG -f $1 -gL $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER

if [ "$4" = "S" ];then   
   echo "/archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap" $DB_37_1 $5".syn " $5".fna"
   time /archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap $DB_37_1 $5.syn $5.fna
fi

echo "Done...!!!"
