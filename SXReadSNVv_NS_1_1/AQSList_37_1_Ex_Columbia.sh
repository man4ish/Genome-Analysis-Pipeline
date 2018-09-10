#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename> <Sample-ID> <Parent folder of raw_read folder>"

if [ $# -ne 3 ]; then
    echo 
    echo $ERRMSG 
    echo $USAGE
    echo "Example:- ./AQSList_37_1.sh ../SRR070490_1_Raw.bin.filtered SRR070490 EXOME_COLUMBIA/SRR070490" 
    echo 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
#READSDIR=/archive/project/oneclick/syhwah/GS000117/NA19240_v36_3/raw_read
READSDIR=/archive/project/oneclick/syhwah/$3/raw_read
REFDIR=/archive/project/oneclick/syhwah/Genome_37.1

FDBSNP=$REFDIR/dbSNP_37.1
FDBINDEL=$REFDIR/dbINDEL_37.1

PRG="/archive/project/oneclick/syhwah/SXReadSNVv_NS_1_1/./SXReadSNV"


echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 30" $2 37.1 33
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 30 $2 37.1 33

echo "Done...!!!"
