#!/bin/sh
export MALLOC_CHECK_=0
READSDIR=/archive/project/oneclick/syhwah/GS000117/37_1_"$3"/raw_read
REFDIR=/archive/project/oneclick/syhwah/Genome_37.1

FDBSNP=$REFDIR/Snp_dbSNP_37.1
FDBINDEL=$REFDIR/Indel_dbSNP_37.1

PRG="/archive/project/oneclick/syhwah/SXReadMultiBpSNV/./SXReadMultiBpSNV"


echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $2 "37.1"
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $2 "37.1"

echo "Done...!!!"
