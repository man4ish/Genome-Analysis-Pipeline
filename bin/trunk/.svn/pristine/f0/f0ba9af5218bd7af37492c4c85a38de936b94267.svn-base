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
    echo "Example:- ./MBSNVList_36_3.sh NS15492_Test.lst 1 500 S NS15492_Test"
    echo  
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
GENOME_VER="36.3"
REFDIR="/archive/project/oneclick/syhwah/Genome_"$GENOME_VER

FBINNED=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.bin
FEXCEPT=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt.exception
FANNOTATE=$REFDIR/1ClickDB_Features_PrimaryAssembly_36.3_250610_Synamatix.txt

PRG="/archive/project/oneclick/syhwah/SXReadSNVIlmn_v1/./SXReadSNVIlmn"

echo $PRG "-f" $1 "-gML" $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER
time $PRG -f $1 -gML $FBINNED $FEXCEPT $FANNOTATE $2 $3 $4 $5 $GENOME_VER

if [ "$4" = "S" ];then
   echo "/archive/project/oneclick/syhwah/SXMBSynoMap/./SXMBSynoMap" $DB_36_3 $5".syn" $5".fna"
   time /archive/project/oneclick/syhwah/SXMBSynoMap/./SXMBSynoMap $DB_36_3 $5.syn $5.fna
fi


echo "Done...!!!"
