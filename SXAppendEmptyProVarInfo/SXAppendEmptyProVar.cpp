#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char** argv)
{

    if (argc !=3){

        fprintf(stdout,"Invalid argument\n");
        fprintf(stdout,"Usage :\t<Input filename> <Output filename>\n");
    }

    FILE *pfIn = fopen(argv[1],"r"); if (!pfIn){fprintf(stdout,"Failed to open %s...\n",argv[1]); exit(9);}
    FILE *pfOut = fopen(argv[2],"w"); if (!pfOut){fprintf(stdout,"Failed to open %s...\n",argv[2]); exit(9);} 
    /*
    fprintf(pfOut,"Chromosome\tGiNumber\tOffset\tNucleotide_Variant\tFwd_SNP_Reads\tRvs_SNP_Reads\t"
                  "Total_SNP_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\t"
                  "Gene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon\tZygosity\t"
                  "Local_Copy_Number\tMappable_Bases\tMappable_Read_Densities\tMappable_Repeat_Densities\t"
                  "Mappable_Bases_Repeat\tRepeat_Densities_MBP\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts)\t"
                  "Synonymous(S)_Non_Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\tPEScore\tSNPReadScore\t"
                  "StrandScore\tReadDenScore\tReptDenScore\tAQSScore\tConfidentScore\tEvidence\tAnnotation\tDrugs\t"
                  "Drug Classes\tDisease\n");
    */

 
    char acbuf[40960], *pChr=NULL;

    while(!feof(pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),pfIn)) break;
     
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';} 
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == '\t'){
           continue;
        }

        if (acbuf[0] == 'C') {
           fprintf(pfOut,"%s\n",acbuf); continue;
        }

        fprintf(pfOut,"%s\t-\t-\t-\t-\n",acbuf); //Chromosome
    } 
   
    if (pfIn) fclose(pfIn); if (pfOut) fclose(pfOut);
}

