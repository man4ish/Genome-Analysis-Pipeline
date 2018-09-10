#!/bin/sh

export MALLOC_CHECK_=0
export M_GENOME="37.1"
export M_CG_SAMPLEID="D01"

APPDIR="/archive/project/oneclick/syhwah"

if [ $M_GENOME = "37.1" ]; then
   READSDIR=$APPDIR"/GS000117/37_1_"$M_CG_SAMPLEID"/raw_read"
elif [ $M_GENOME = "36.3" ]; then
   READSDIR=$APPDIR"/GS000117/"$M_CG_SAMPLEID"/run2/raw_read"
else
   READSDIR=$APPDIR"/Illumina/Harvard/run_2/raw_read"
fi

REFDIR=$APPDIR"/Genome_"$M_GENOME
FDBSNP=$REFDIR"/Snp_dbSNP_"$M_GENOME
FDBINDEL=$REFDIR"/Indel_dbSNP_"$M_GENOME

PRG=$APPDIR"/SXReadSNVv4/./SXReadSNV"

echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $2
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $2

export M_SAMPLEID=$2
export M_FILTER_LIST=""     #value will be defined in ProcSnp.sh
export M_MIN_PSR="1"        #min supporting reads to be printed in final stats report
export M_MAX_PSR="500"      #max supporting reads to be printed in final stats report
export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del  


#Process SNP List ...
. $APPDIR"/SXReadSNVv4/FilterSNV.sh"
. $APPDIR"/SXReadSNVv4/GMasterList.sh"

M_TYPE="I" 

#Process INS List ...
. $APPDIR"/SXReadSNVv4/FilterSNV.sh"
. $APPDIR"/SXReadSNVv4/GMasterList.sh"

M_TYPE="D"

#Process DEL List ...
. $APPDIR"/SXReadSNVv4/FilterSNV.sh"
. $APPDIR"/SXReadSNVv4/GMasterList.sh"


echo "Done...!!!"
