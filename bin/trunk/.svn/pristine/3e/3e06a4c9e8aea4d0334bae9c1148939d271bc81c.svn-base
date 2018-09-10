#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename> <Sample-ID> <Parent folder of raw_read folder>"

if [ $# -ne 3 ]; then
    echo 
    echo $ERRMSG 
    echo $USAGE
    echo "Example:- ./MBAQSList_36_3.sh ../GS000117/NS15492_TEST/TEST_21_gbl_mnv_abs NS15492_TEST_MB NS15492_TEST" 
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

PRG="/archive/project/oneclick/syhwah/SXReadSNVIlmn_v1/./SXReadSNVIlmn"


echo $PRG "-f" $1 "-gMLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 30" $2 36.3 64
time $PRG -f $1 -gMLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 30 $2 36.3 64

echo "Done...!!!"
