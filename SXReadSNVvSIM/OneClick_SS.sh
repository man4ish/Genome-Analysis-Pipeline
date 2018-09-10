#!/bin/sh

#-----------------------------------------------------------------------------------------------
ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename> <Sample-ID> <Genome Ver.>" #<Mean*3sd> <Mean Value 1> <Min. INDEL Zygo SR>"

if [ $# -ne 3 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------


export MALLOC_CHECK_=0
export M_GENOME="36.3"
export M_CG_SAMPLEID="NA19240"
export M_APPDIR="/archive/project/oneclick/syhwah"

if [ $M_GENOME = "37.1" ]; then
   READSDIR=$M_APPDIR"/GS000117/37_1_"$M_CG_SAMPLEID"/raw_read"
elif [ $M_GENOME = "36.3" ]; then
   READSDIR=$M_APPDIR"/GS000117/"$M_CG_SAMPLEID"/raw_read"
else
   READSDIR=$M_APPDIR"/Illumina/Harvard/run_2/raw_read"
fi

REFDIR=$M_APPDIR"/Genome_"$M_GENOME
FDBSNP=$REFDIR"/dbSNP_"$M_GENOME
FDBINDEL=$REFDIR"/dbINDEL_"$M_GENOME
PRG=$M_APPDIR"/SXReadSNVv8_1_2/./SXReadSNV"

echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINDEL $READSDIR "1 10000 0 0" $2 $M_GENOME 64
#time $PRG -f $1 -gLa $FDBSNP $FDBINDEL $READSDIR 1 10000 0 0 $2 $M_GENOME 64

export M_MIN_KSR=2           #Min. Known Supporting Reads
export M_MIN_NSR=4           #Min. Novel Supporting Reads               
export M_MIN_NRS=15          #Min. Read Strength
export M_MIN_NRD=4           #Min. Read Density
export M_MAX_SR=180          #Max. Known/Novel Supporting Reads

export M_SAMPLEID=$2
#export M_Mean3SD=$3
#export M_Mean=$4
#export M_Zygo_MinSR=$5

export M_FILTER_LIST=""     #value will be defined in ProcSnp.sh
export M_MIN_PSR="1"        #min supporting reads to be printed in final stats report
export M_MAX_PSR="500"      #max supporting reads to be printed in final stats report
export M_ZYGO_CUTOFF="74"     #Zygosity cut-off point
export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del  


#Process SNP List ...
. $M_APPDIR"/SXReadSNVv8_1_2/OneClickFilter_SS.sh"

#export M_SNP_ZYGO_CUTOFF=72
#export M_MEAN=$M_Mean
#export M_SNP_HET_MSR=4
#export M_SNP_HOM_MSR=5

. $M_APPDIR"/SXReadSNVv8_1_2/OneClickAnnotate_SS.sh"


#Process INS List ...

M_TYPE="I"
M_MIN_KSR=1
M_MIN_NSR=3
M_MIN_NRS=15
M_MIN_NRD=3

. $M_APPDIR"/SXReadSNVv8_1_2/OneClickFilter_SS.sh"
. $M_APPDIR"/SXReadSNVv8_1_2/OneClickAnnotate_SS.sh"


#Process DEL List ...

M_TYPE="D"

. $M_APPDIR"/SXReadSNVv8_1_2/OneClickFilter_SS.sh"
. $M_APPDIR"/SXReadSNVv8_1_2/OneClickAnnotate_SS.sh"


