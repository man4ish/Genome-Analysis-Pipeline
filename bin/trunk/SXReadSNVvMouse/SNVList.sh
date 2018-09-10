#!/bin/sh
export MALLOC_CHECK_=0

REFDIR=/archive/project/oneclick/syhwah/Genome_37.1

FGENE=$REFDIR/Genedb_37.1
FEXON=$REFDIR/exon_37.1_ccds #CCDS_Exon
FCNV=$REFDIR/cnv_Variation_37.1_ncbi

PRG="/archive/project/oneclick/syhwah/SXReadSNVv4/./SXReadSNV"

echo $PRG "-f" $1 "-gL" $FGENE $FEXON $FCNV $2 $3 $4 $5
time $PRG -f $1 -gL $FGENE $FEXON $FCNV $2 $3 $4 $5

if [ "$4" = "S" ];then
   DB=$DB_37_1 

   echo "/archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap" $DB $5".syn " $5".fna"
   time /archive/project/oneclick/syhwah/SXSynoMap/./SXSynoMap $DB $5.syn $5.fna
fi

echo "Done...!!!"
