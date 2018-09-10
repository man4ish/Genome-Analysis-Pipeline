#!/bin/sh
export MALLOC_CHECK_=0
#READSDIR=/archive/project/oneclick/syhwah/GS000117/NA19240_v36_3/raw_read
READSDIR=/archive/project/oneclick/syhwah/GS000117/$3/raw_read
REFDIR=/archive/project/oneclick/syhwah/Genome_36.3

FDBSNP=$REFDIR/dbSNP_36.3
FDBINDEL=$REFDIR/dbINDEL_36.3

PRG="/archive/project/oneclick/syhwah/SXReadSNVv8_1_1/./SXReadSNV"


echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $2
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $2

echo "Done...!!!"
