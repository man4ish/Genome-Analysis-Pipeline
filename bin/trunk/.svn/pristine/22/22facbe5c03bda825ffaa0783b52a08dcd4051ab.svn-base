#!/bin/sh
export MALLOC_CHECK_=0
READSDIR=/archive/project/oneclick/syhwah/ANU1006/raw_read
REFDIR=/archive/project/oneclick/syhwah/Mouse_Genome

FDBSNP=$REFDIR/dbSNP128_Mouse_SNP
FDBINDEL=$REFDIR/dbSNP128_Mouse_INDEL

PRG="/archive/project/oneclick/syhwah/SXReadSNVvMouse/./SXReadSNV"

echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 25" $2
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 25 $2

echo "Done...!!!"
