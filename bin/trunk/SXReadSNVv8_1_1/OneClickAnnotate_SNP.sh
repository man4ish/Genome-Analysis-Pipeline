#!/bin/sh

#-----------------------------------------------------------------------------------------------

ERRMSG="Invalid Parameters..."
USAGE="Usage: `basename $0` <SNV Input filename 1> <Zygosity CutOff> <Mean Value> <Sample-ID>"

if [ $# -ne 4 ]; then
    echo $ERRMSG 
    echo $USAGE
    echo "" 
    exit 1
fi

#------------------------------------------------------------------------------------------------

FILE_INPUT=$1
ZYGO_CUTOFF=$2
MEANS=$3
TYPE="S"
SAMPLEID=$4

HET_MSR=4
HOM_MSR=5


M_GENOME="37.1"
#M_CG_SAMPLEID="NA19240"
M_APPDIR="/archive/project/oneclick/syhwah"

READSDIR=$M_APPDIR"/GS000117/37_1_"$SAMPLEID"/raw_read"
REPTDENDIR=$READSDIR
MASKDENDIR=$M_APPDIR"/SXFilterRepeat/Reptmask"

REFDIR=$M_APPDIR"/Genome_"$M_GENOME

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV" 

echo $PRG $FILE_INPUT $ZYGO_CUTOFF $HET_MSR $HOM_MSR $MEANS $READSDIR $REPTDENDIR $MASKDENDIR $FILE_INPUT"_ZygoCNV"
time $PRG $FILE_INPUT $ZYGO_CUTOFF $HET_MSR $HOM_MSR $MEANS $READSDIR $REPTDENDIR $MASKDENDIR $FILE_INPUT"_ZygoCNV"

file_size=(`du -cb $FILE_INPUT"_ZygoCNV"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FBINNED=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVv8_1_1/./SXReadSNV"

echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $TYPE $SAMPLEID $M_GENOME
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $TYPE $SAMPLEID $M_GENOME

if [ "$M_GENOME" = "37.1" ]; then  
   DB=$DB_37_1 
else
   DB=$DB_36_3
fi

echo $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"
time $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $SAMPLEID".syn" $SAMPLEID".fna"

file_size=(`du -cb $SAMPLEID".fna"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#echo "Done...!!!"

#------------------------------------------------------------------------------------------
#synasearch 

#echo $SYNASEARCH -np -f $SAMPLEID".fna" -db $DB > $SAMPLEID".sso"
#time eval $SYNASEARCH -np -f $SAMPLEID".fna" -db $DB > $SAMPLEID".sso"
#
#file_size=(`du -cb $SAMPLEID".sso"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------

CURR_DATE=$(date '+%y%m%d')
PRG=$M_APPDIR"/SSNpParser/./SSNpParser" 
IN_FILE=$SAMPLEID"_snp_c1_"$CURR_DATE
OUT_FILE=$IN_FILE".lst"

echo $PRG $SAMPLEID".sso" $IN_FILE $OUT_FILE $SAMPLEID  
time $PRG $SAMPLEID".sso" $IN_FILE $OUT_FILE $SAMPLEID

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendScore/./SXAppendScore"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".score"

echo $PRG $IN_FILE $OUT_FILE
time $PRG $IN_FILE $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendDrugInfo/./SXAppendDrugInfo"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".druginfo"

echo $PRG "dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE
time $PRG $M_APPDIR"/Genome_37.1/dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi



