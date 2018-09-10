#!/bin/sh

#-----------------------------------------------------------------------------------------------

ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename 1> <Sample-ID 1> <SNV Input filename 2> \
                            <Sample-ID 2> <Mean*3sd Value 1> <Mean*3sd Value 2> \
                            <Avg. Mean*3sd of Value 1 Value 2> <Mean Valiue 1> \
                            <Mean Value 2> <Min. INDEL Zygo SR 1>  <Min. INDEL Zygo SR 2>"

if [ $# -ne 11 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export M_GENOME="37.1"
export M_APPDIR="/archive/project/oneclick/syhwah"
export M_MASKDENSDIR=$M_APPDIR"/SXFilterRepeat/Reptmask"

export M_REFDIR=$M_APPDIR"/Genome_"$M_GENOME
FDBSNP=$M_REFDIR"/dbSNP_131" #$M_GENOME
FDBINS=$M_REFDIR"/dbINS_131" #$M_GENOME
FDBDEL=$M_REFDIR"/dbDEL_131" #$M_GENOME
HSREF="hs_ref_37.1"
#------------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXReadSNVvPoly2_dbSNP131/./SXReadSNV"

READSDIR=$M_APPDIR"/GS000117/37_1_"$2"/raw_read"
echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINS $FDBDEL $READSDIR "1 10000 0 0" $2 $M_GENOME "33" $HSREF
#time $PRG -f $1 -gLa $FDBSNP $FDBINS $FDBDEL $READSDIR 1 10000 0 0 $2 $M_GENOME "33" $HSREF


READSDIR=$M_APPDIR"/GS000117/37_1_"$4"/raw_read"
echo $PRG "-f" $3 "-gLa" $FDBSNP $FDBINS $FDBDEL $READSDIR "1 10000 0 0" $4 $M_GENOME 33 $HSREF
#time $PRG -f $3 -gLa $FDBSNP $FDBIN $FDBDEL $READSDIR 1 10000 0 0 $4 $M_GENOME 33 $HSREF

#------------------------------------------------------------------------------------------------

export M_MIN_KSR=3           #Min. Known Supporting Reads
export M_MIN_NSR=3           #Min. Novel Supporting Reads               
export M_MIN_NRS=20          #Min. Read Strength
export M_MIN_NRD=3           #Min. Read Density

export M_MIN_COMMON_KSR=3
export M_MIN_COMMON_NSR=3
export M_MIN_COMMON_NRD=3

export M_MIN_KSR2=3           #Min. Known Supporting Reads
export M_MIN_NSR2=3           #Min. Novel Supporting Reads               
export M_MIN_NRD2=3

export M_MAX_SR=180          #Max. Known/Novel Supporting Reads

export M_SAMPLEID_1=$2
export M_SAMPLEID_2=$4

export M_Mean3SD_1=$5
export M_Mean3SD_2=$6
export M_Mean3SD_Common=$7
export M_Mean_1=$8
export M_Mean_2=$9

shift 2 
export M_Zygo_MinSR_1=$8
export M_Zygo_MinSR_2=$9

export M_TYPE="S"           #type=> S-Snp, I-Ins, D-Del 

#------------------------------------------------------------------------------------------------

#Process SNP List ...
   
#. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickFilter_IS.sh"

export M_SAMPLEID=$M_SAMPLEID_1    
export M_SNP_ZYGO_CUTOFF=72
export M_MEAN=$M_Mean_1
export M_SNP_HET_MSR=4
export M_SNP_HOM_MSR=5

#. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS_SNP.sh"


M_SAMPLEID=$M_SAMPLEID_2
M_MEAN=$M_Mean_2

#. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS_SNP.sh"
#------------------------------------------------------------------------------------------------

#Process INS List ...

M_TYPE="I"

M_MIN_KSR=6           #Min. Known Supporting Reads
M_MIN_NSR=6           #Min. Novel Supporting Reads               
M_MIN_NRD=6

M_MIN_KSR2=6           #Min. Known Supporting Reads
M_MIN_NSR2=6           #Min. Novel Supporting Reads               
M_MIN_NRD2=6

M_MIN_COMMON_KSR=6
M_MIN_COMMON_NSR=6
M_MIN_COMMON_NRD=6

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickFilter_IS.sh"


M_SAMPLEID=$M_SAMPLEID_1
export M_SAMPLENO=1
export M_DIPLOID_PROB=5
M_MEAN=$M_Mean_1
export M_ZYGO_MIN_SR=$M_Zygo_MinSR_1

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS.sh"

M_SAMPLEID=$M_SAMPLEID_2
M_SAMPLENO=2
M_MEAN=$M_Mean_2
M_ZYGO_MIN_SR=$M_Zygo_MinSR_2

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS.sh"

#------------------------------------------------------------------------------------------------

#Process DEL List ...

M_TYPE="D"

M_MIN_KSR=6           #Min. Known Supporting Reads
M_MIN_NSR=6           #Min. Novel Supporting Reads               
M_MIN_NRD=6 

M_MIN_KSR2=6           #Min. Known Supporting Reads
M_MIN_NSR2=6           #Min. Novel Supporting Reads               
M_MIN_NRD2=6

M_MIN_COMMON_KSR=6
M_MIN_COMMON_NSR=6
M_MIN_COMMON_NRD=6

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickFilter_IS.sh"

M_SAMPLEID=$M_SAMPLEID_1
M_SAMPLENO=1
M_MEAN=$M_Mean_1
M_ZYGO_MIN_SR=$M_Zygo_MinSR_1

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS.sh"

M_SAMPLEID=$M_SAMPLEID_2
M_SAMPLENO=2
M_MEAN=$M_Mean_2
M_ZYGO_MIN_SR=$M_Zygo_MinSR_2

. $M_APPDIR"/SXReadSNVvPoly2_dbSNP131/OneClickAnnotate_IS.sh"

