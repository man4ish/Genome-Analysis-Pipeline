#!/bin/sh

#-----------------------------------------------------------------------------------------------

FILE_INPUT=$M_SAMPLEID"_snp_Sample"

READSDIR=$M_APPDIR"/GS000117/37_1_"$M_SAMPLEID"/raw_read"
REPTDENDIR=$READSDIR

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintZygosityCNV/./SXPrintZygosityCNVEx" 
                                                                                                       
echo $PRG $FILE_INPUT $M_SNP_ZYGO_CUTOFF $M_SNP_HET_MSR $M_SNP_HOM_MSR $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $M_SAMPLEID $FILE_INPUT"_ZygoCNV"
time $PRG $FILE_INPUT $M_SNP_ZYGO_CUTOFF $M_SNP_HET_MSR $M_SNP_HOM_MSR $M_MEAN $READSDIR $REPTDENDIR $M_MASKDENSDIR $M_SAMPLEID $FILE_INPUT"_ZygoCNV"

file_size=(`du -cb $FILE_INPUT"_ZygoCNV"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

FBINNED=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.bin
FEXCEPT=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt.exception
FANNOTATE=$M_REFDIR/GRCh37.1_H-Sapiens-010710-1ClickDB_Features-Synamatix.txt

PRG=$M_APPDIR"/SXReadSNV_CG/./SXReadSNV_CG"

echo $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME $M_OUTPUT_DIR
time $PRG "-f" $FILE_INPUT"_ZygoCNV" "-gL" $FBINNED $FEXCEPT $FANNOTATE 1 500 $M_TYPE $M_SAMPLEID $M_GENOME $M_OUTPUT_DIR


#------------------------------------------------------------------------------------------

if [ -z "$SYNABASEROOT" ] ; then
  if [ -d "/opt/synamatix/synabase/synabase" ]; then
    export SYNABASEROOT=/opt/synamatix/synabase/synabase
    export LD_LIBRARY_PATH=$SYNABASEROOT/lib:$LD_LIBRARY_PATH
  else
    echo "The SYNABASEROOT environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SXDBConfig" ] ; then
  if [ -f "/data/dbref/synabase/sxdbconfig.xml" ]; then
    export SXDBConfig=/archive/dbref/synabase/sxdbconfig.xml
  else
    echo "The SXDBConfig environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SYNASEARCH" ] ; then
  if [ -f "/opt/synamatix/synasearch/synasearch/bin/synasearch.sh" ]; then
    export SYNASEARCH=/opt/synamatix/synasearch/synasearch/bin/synasearch.sh
  else
    echo "The SYNASEARCH environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SXDB_CG" ] ; then
 SXDB_CG="sprt_prtn_upkb-57.8_sb3.0.9"
fi
#------------------------------------------------------------------------------------------


echo $M_APPDIR"/SXSynoMap/./SXSynoMapEx" $M_HSREF $M_GENOME $M_SAMPLEID".syn" $M_SAMPLEID".fna"
time $M_APPDIR"/SXSynoMap/./SXSynoMapEx" $M_HSREF $M_GENOME $M_SAMPLEID".syn" $M_SAMPLEID".fna"

file_size=(`du -cb $M_SAMPLEID".fna"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
#synasearch 

echo $SYNASEARCH "-np -f" $M_SAMPLEID".fna -db" $SXDB_CG ">" $M_SAMPLEID".sso"
time $SYNASEARCH -np -f $M_SAMPLEID".fna" -db $SXDB_CG > $M_SAMPLEID".sso"

file_size=(`du -cb $M_SAMPLEID".sso"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
#
CURR_DATE=$(date '+%y%m%d')
PRG=$M_APPDIR"/SSNpParser/./SSNpParser" 
IN_FILE=$M_SAMPLEID"_snp_c1_"$CURR_DATE
OUT_FILE=$IN_FILE".lst"

echo $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID  
time $PRG $M_SAMPLEID".sso" $IN_FILE $OUT_FILE $M_SAMPLEID

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendScore/./SXAppendScoreSNP"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".score"

echo $PRG $IN_FILE $M_MEAN $OUT_FILE
time $PRG $IN_FILE $M_MEAN $OUT_FILE

file_size=(`du -cb $OUT_FILE`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXPrintStats/SXPrintStatsSNP"

echo  $PRG $OUT_FILE "1" $OUT_FILE".stats.rpt" 
time  $PRG $OUT_FILE "1" $OUT_FILE".stats.rpt" 

file_size=(`du -cb $OUT_FILE".stats.rpt"`)
if [ $file_size -eq 0 ]; then
   exit 1
fi

#------------------------------------------------------------------------------------------
#Update Summary Report

PRG=$M_APPDIR"/Avg_Confi_Score/Avg_Confi_Score"

SUMMARY_RPT=$M_SAMPLEID"_snp_summary_c1_"$CURR_DATE".rpt"
SUMMARY_ORI=$M_SAMPLEID"_snp_summary_c1_"$CURR_DATE".rpt.ori"

mv $SUMMARY_RPT $SUMMARY_ORI

echo $PRG $OUT_FILE $SUMMARY_ORI $SUMMARY_RPT
time $PRG $OUT_FILE $SUMMARY_ORI $SUMMARY_RPT

#------------------------------------------------------------------------------------------

PRG=$M_APPDIR"/SXAppendDrugInfo/./SXAppendDrugInfo"
IN_FILE=$OUT_FILE
OUT_FILE=$IN_FILE".druginfo"

echo $PRG "dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE
time $PRG $M_APPDIR"/Genome_37.1/dbSNP131_37.1_hg19_17August_Synamatix.txt" $IN_FILE $OUT_FILE

#file_size=(`du -cb $OUT_FILE`)
#if [ $file_size -eq 0 ]; then
#   exit 1
#fi



