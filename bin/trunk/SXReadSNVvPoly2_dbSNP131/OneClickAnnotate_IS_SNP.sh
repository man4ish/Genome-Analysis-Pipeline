#!/bin/sh

#-----------------------------------------------------------------------------------------------

FILE_INPUT=$M_SAMPLEID"_snp_Sample"

READSDIR=$M_APPDIR"/GS000117/37_1_"$M_SAMPLEID"/raw_read"
REPTDENDIR=$READSDIR

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNV" 
                                                                                                       
echo $PRG $FILE_INPUT $M_SNP_ZYGO_CUTOFF $M_SNP_HET_MSR $M_SNP_HOM_MSR $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $FILE_INPUT"_ZygoCNV"
time $PRG $FILE_INPUT $M_SNP_ZYGO_CUTOFF $M_SNP_HET_MSR $M_SNP_HOM_MSR $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $FILE_INPUT"_ZygoCNV"

file_size=(`du -cb $FILE_INPUT"_ZygoCNV"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FBINNED=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNVvPoly2_dbSNP131/./SXReadSNV"

echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME

if [ "$M_GENOME" = "37.1" ]; then  
   DB=$DB_37_1 
else
   DB=$DB_36_3
fi

echo $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $M_SAMPLEID".syn" $M_SAMPLEID".fna"
time $M_APPDIR"/SXSynoMap/./SXSynoMap" $DB $M_SAMPLEID".syn" $M_SAMPLEID".fna"

file_size=(`du -cb $M_SAMPLEID".fna"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#echo "Done...!!!"

#------------------------------------------------------------------------------------------
#synasearch 

#echo $SYNASEARCH -np -f $M_SAMPLEID".fna" -db $DB > $M_SAMPLEID".sso"
#time eval $SYNASEARCH -np -f $M_SAMPLEID".fna" -db $DB > $M_SAMPLEID".sso"
#
#file_size=(`du -cb $M_SAMPLEID".sso"`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------
#
#CURR_DATE=$(date '+%y%m%d')
#PRG=$M_APPDIR"/SSNpParser/./SSNpParser" 
#IN_FILE=$M_SAMPLEID"_snp_c1_"$CURR_DATE
#OUT_FILE=$IN_FILE".lst"
#
#echo $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID  
#time $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID
#
#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreSNP"
#IN_FILE=$OUT_FILE
#OUT_FILE=$IN_FILE".score"

#echo $PRG $IN_FILE $M_MEAN $OUT_FILE
#time $PRG $IN_FILE $M_MEAN $OUT_FILE

#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi

#------------------------------------------------------------------------------------------

#PRG=$M_APPDIR"/SXAppendDrugInfo/./SXAppendDrugInfo"
#IN_FILE=$OUT_FILE
#OUT_FILE=$IN_FILE".druginfo"

#echo $PRG "dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE
#time $PRG $M_APPDIR"/Genome_37.1/dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE

#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi



