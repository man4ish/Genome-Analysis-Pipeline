#!/bin/sh

#------------------------------------------------------------------------------------------------

exporRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <Sample-ID> <Raw Reads folder Path>" 

if [ $# -ne 11 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------i
export M_GENOME="37.1"
export M_CG_SAMPLEID=$1
export M_APPDIR="/archive/project/oneclick/syhwah"

export M_READSDIR=$M_APPDIR"/"$2"/raw_read"

REFDIR=$M_APPDIR"/Genome_"$M_GENOME
FDBSNP=$REFDIR"/dbSNP_"$M_GENOME
FDBINDEL=$REFDIR"/dbINDEL_"$M_GENOME
#PRG=$M_APPDIR"/SXReadSNVv6/./SXReadSNV"

#echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $M_READSDIR "1 10000 0 0" $2
#time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $M_READSDIR 1 10000 0 0 $2

export M_MIN_KSR=3           #Min. Known Supporting Reads
export M_MIN_NSR=3           #Min. Novel Supporting Reads               
export M_MIN_NRS=15          #Min. Read Strength
export M_MIN_NRD=3           #Min. Read Density
export M_MAX_SR=100          #Max. Known/Novel Supporting Reads

export M_SAMPLEID=$1
export M_FILTER_LIST=""     #value will be defined in ProcSnp.sh
export M_MIN_PSR="1"        #min supporting reads to be printed in final stats report
export M_MAX_PSR="500"      #max supporting reads to be printed in final stats report
export M_ZYGO_CUTOFF="72"    #Zygosity cut-off point
export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del  


#Process SNP List ...
. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickFilter_SS.sh"
. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickAnnotate_SS_SNP_EXOME_DEMO.sh"


#Process INS List ...

M_TYPE="I"
M_MIN_KSR=1
M_MIN_NSR=3
M_MIN_NRS=15
M_MIN_NRD=3
M_ZYGO_MIN_SR=3
M_DIPLOID_PROB=5
M_SAMPLENO=1

. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickFilter_SS.sh"
. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickAnnotate_SS_EXOME_DEMO.sh"


#Process DEL List ...

M_TYPE="D"

. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickFilter_SS.sh"
. $M_APPDIR"/SXReadSNVv_NS_1_1/OneClickAnnotate_SS_EXOME_DEMO.sh"


