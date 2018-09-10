#!/bin/sh

#-----------------------------------------------------------------------------------------------

ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename 1> <Sample-ID 1> <SNV Input filename 2> \
                            <Sample-ID 2> <Mean*3sd Value 1> <Mean*3sd Value 2> \
                            <Avg. Mean*3sd of Value 1 Value 2>"

if [ $# -ne 7 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export M_GENOME="37.1"
export M_APPDIR="/archive/project/oneclick/syhwah"
export M_MASKDENSDIR=$M_APPDIR"/SXFilterRepeat/Reptmask"

REFDIR=$M_APPDIR"/Genome_"$M_GENOME
FDBSNP=$REFDIR"/dbSNP_"$M_GENOME
FDBINDEL=$REFDIR"/dbINDEL_"$M_GENOME

#------------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXReadSNVv8_1_1/./SXReadSNV"

READSDIR=$M_APPDIR"/GS000117/37_1_"$2"/raw_read"
echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $2
time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $2


READSDIR=$M_APPDIR"/GS000117/37_1_"$4"/raw_read"
echo $PRG "-f" $3 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $4
time $PRG -f $3 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $4

#------------------------------------------------------------------------------------------------

export M_MIN_KSR=2           #Min. Known Supporting Reads
export M_MIN_NSR=3           #Min. Novel Supporting Reads               
export M_MIN_NRS=20          #Min. Read Strength
export M_MIN_NRD=3           #Min. Read Density
export M_MAX_SR=180          #Max. Known/Novel Supporting Reads

export M_SAMPLEID_1=$2
export M_SAMPLEID_2=$4

export M_Mean3SD_1=$5
export M_Mean3SD_2=$6
export M_Mean3SD_Common=$7
export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del 

#------------------------------------------------------------------------------------------------

#Process SNP List ...
   

. $M_APPDIR"/SXReadSNVv8_1_1/OneClickFilter.sh"
#. $M_APPDIR"/SXReadSNVv6/OneClickFinal.sh"

#------------------------------------------------------------------------------------------------

#Process INS List ...

M_TYPE="I"

. $M_APPDIR"/SXReadSNVv8_1_1/OneClickFilter.sh"
#. $M_APPDIR"/SXReadSNVv6/OneClickFinal.sh"

#------------------------------------------------------------------------------------------------

#Process DEL List ...

M_TYPE="D"

. $M_APPDIR"/SXReadSNVv8_1_1/OneClickFilter.sh"
#. $M_APPDIR"/SXReadSNVv6/OneClickFinal.sh"


