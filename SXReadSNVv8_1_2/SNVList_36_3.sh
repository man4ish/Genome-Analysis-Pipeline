#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <Input filename> <Min Supporting Reads to be printed> \
                            <Max Supporting Reads to be printed> \
                            <SNV Type eg. S,I,D> <Sample-ID>"

if [ $# -ne 5 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
GENOME_VER="36.3"
REFDIR="/archive/project/oneclick/syhwah/Genome_"$GENOME_VER

FBINNED=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
FEXCEPT=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
FANNOTATE=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt

PRG="/archive/project/oneclick/syhwah/SXReadSNVv8_1_2/./SXReadSNV"

echo $PRG "-f" $1 "-gL" $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER
time $PRG -f $1 -gL $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER

if [ "$4" = "S" ];then
   echo "/archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap" $DB_36_3 $5".syn " $5".fna"
   time /archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap $DB_36_3 $5.syn $5.fna
fi


echo "Done...!!!"
