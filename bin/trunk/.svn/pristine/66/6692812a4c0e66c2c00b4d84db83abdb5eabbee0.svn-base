#!/bin/sh

#-----------------------------------------------------------------------------------------------

ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename 1> <Sample-ID 1> <Ambiguity filename 1>" 

if [ $# -ne 3 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

export M_GENOME="37.1"
export M_APPDIR="/archive/project/oneclick/syhwah"

export M_REFDIR=$M_APPDIR"/Genome_"$M_GENOME
FDBSNP=$M_REFDIR"/dbSNP_131"
FDBINS=$M_REFDIR"/dbINS_131"
FDBDEL=$M_REFDIR"/dbDEL_131"
export M_HSREF="../hs_ref_37.1"
M_OUTPUT_DIR="."
PROJ_SAMPLEID=$2
#------------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXReadSNV_CG/./SXReadSNV_CG"

READSDIR=$M_APPDIR"/GS000117/37_1_"$2"/raw_read"
echo $PRG "-f" $1 "-gLa" $FDBSNP $FDBINS $FDBDEL $READSDIR "1 100000 0 0" $2 $M_GENOME "33" \
     $M_HSREF $M_OUTPUT_DIR $PROJ_SAMPLEID "-a" $3
time $PRG -f $1 -gLa $FDBSNP $FDBINS $FDBDEL $READSDIR 1 100000 0 0 $2 $M_GENOME 33 $M_HSREF \
     $M_OUTPUT_DIR $PROJ_SAMPLEID -a $3

