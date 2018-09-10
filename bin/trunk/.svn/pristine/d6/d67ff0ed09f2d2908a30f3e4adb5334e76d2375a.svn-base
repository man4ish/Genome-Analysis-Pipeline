#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename> <Sample-ID> <Parent folder of raw_read folder> 
       <Human Ref. files path> <SNV Ambiguity filename>"

if [ $# -ne 5 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export MALLOC_CHECK_=0
READSDIR=/archive/project/oneclick/syhwah/GS000117/$3/raw_read
REFDIR=/archive/project/oneclick/syhwah/Genome_37.1

FDBSNP=$REFDIR/Snp_dbSNP_37.1
FDBINDEL=$REFDIR/Indel_dbSNP_37.1

PRG="/archive/project/oneclick/syhwah/SXReadSNVvSIMPolyPDens2/./SXReadSNV"


echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 100000 0 0" $2 37.1 33 $4 $5
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 100000 0 0 $2 37.1 33 $4 $5
echo "Done...!!!"
