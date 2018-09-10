#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename> <Sample-ID> <Parent folder of raw_read folder>"

if [ $# -ne 3 ]; then
    echo 
    echo $ERRMSG 
    echo $USAGE
    echo "Example:- ./AQSList_36_3.sh ../NS15492_1_raw.bin.filtered NS15492 NS15492" 
    echo 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
#READSDIR=/archive/project/oneclick/syhwah/GS000117/NA19240_v36_3/raw_read
READSDIR=/archive/project/oneclick/syhwah/GS000117/$3/raw_read
REFDIR=/archive/project/oneclick/syhwah/Genome_36.3

FDBSNP=$REFDIR/dbSNP_36.3
FDBINDEL=$REFDIR/dbINDEL_36.3

PRG="/archive/project/oneclick/syhwah/SXReadSNVv_NS_1_1/./SXReadSNV"


echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 30" $2 36.3 64
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 30 $2 36.3 64

echo "Done...!!!"
