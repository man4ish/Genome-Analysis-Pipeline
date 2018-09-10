#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <Input filename> <Min Supporting Reads to be printed> \
                            <Max Supporting Reads to be printed> \
                            <SNV Type eg. S,I,D> <Sample-ID>"

if [ $# -ne 5 ]; then
    echo 
    echo $ERRMSG 
    echo $USAGE
    echo "Example:- ./MBSNVList_37_1.sh MBNS15492_Test_37_1.lst 1 500 S MBNS15492_Test_37_1"
    echo  
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
GENOME_VER="37.1"
REFDIR="/archive/project/oneclick/syhwah/Genome_"$GENOME_VER

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG="/archive/project/oneclick/syhwah/SXReadSNVIlmn_v1/./SXReadSNVIlmn"

echo $PRG "-f" $1 "-gML" $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER
time $PRG -f $1 -gML $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER

if [ "$4" = "S" ];then
   echo "/archive/project/oneclick/syhwah/SXMBSynoMap/./SXMBSynoMap" $DB_37_1 $5".syn" $5".fna"
   time /archive/project/oneclick/syhwah/SXMBSynoMap/./SXMBSynoMap $DB_37_1 $5.syn $5.fna
fi


echo "Done...!!!"
