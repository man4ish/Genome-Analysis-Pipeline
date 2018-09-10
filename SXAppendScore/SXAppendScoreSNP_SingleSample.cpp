#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL, *g_pfOut=NULL; float g_fMean=0.00;


static void banner(char *argv[])
{   
    printf("Synamatix SXAppendScore Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Mean> <Output filename>\n\n");
}

void ProcessRecs()
{
    int nFwdSNPCnt,nRvsSNPCnt,nSNPCnt,nPECnt,nTotalRD,nMappableBases,nMappableReptDens;
    int nReptDenScore,nMin;  
    float fPEScore,fSNPReadScore,fStrandScore,fAvgQS,fRS,fReadDenScore,fAQSScore,fSNPScore;  
    float fXCoverage=g_fMean/2.00, fYCoverage=g_fMean;

    char acbuf[40960],acOut[40960],*pChr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {        
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %lu x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }
    
        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        if (acbuf[0] == 'C') {
           fprintf(g_pfOut,"%s\tPEScore\tSNPReadScore\tStrandScore\tReadDenScore"
                           "\tReptDenScore\tAQSScore\tConfidentScore\n",acbuf);          
           continue;
        }
         
        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t"); pChr = strtok(NULL,"\t"); //Chromosome,GiNumber
        pChr = strtok(NULL,"\t");  pChr = strtok(NULL,"\t"); //Offset,Nucleotide_Variant

        pChr = strtok(NULL,"\t");  nFwdSNPCnt = atoi(pChr); //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  nRvsSNPCnt = atoi(pChr); //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  nSNPCnt = atoi(pChr);    //SNP_Reads                 
        pChr = strtok(NULL,"\t");  nTotalRD = atoi(pChr);   //Total_Read_Density
        pChr = strtok(NULL,"\t");  nPECnt = atoi(pChr);     //PEnd_Count 

        pChr = strtok(NULL,"\t"); //Avg_QryPos 
        pChr = strtok(NULL,"\t"); fAvgQS = atof(pChr);//Avg_QScore
        pChr = strtok(NULL,"\t"); //dbSNP      

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");//Gene_Name,Gene_Description,Keyword  
        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");//miRNA,Promoter,UTR
        pChr = strtok(NULL,"\t"); //Exon 
        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t"); //Zygosity,Local_Copy_Number
        
        pChr = strtok(NULL,"\t"); nMappableBases = atoi(pChr);//Mappable_Bases      

        pChr = strtok(NULL,"\t"); // Mappable_Read_Densities
        pChr = strtok(NULL,"\t"); nMappableReptDens = atoi(pChr);//Mappable_Repeat_Densities

        fPEScore = !(nPECnt < fXCoverage)?100.00:(float)nPECnt*(100.00/fXCoverage);
        fSNPReadScore = !(nSNPCnt < fYCoverage)?100.00:(float)nSNPCnt*(100.00/fYCoverage);

        if (!(nFwdSNPCnt < fXCoverage) && !(nRvsSNPCnt < fXCoverage)) fStrandScore = 100.00;
        else if ((nFwdSNPCnt < 2) && (nRvsSNPCnt < 2)) fStrandScore = 0.00; 
        else {
           nMin = (nFwdSNPCnt < nRvsSNPCnt)? nFwdSNPCnt:nRvsSNPCnt;
           fStrandScore = (float)nMin*(100.00/fXCoverage);                    
        } 

        fRS = (float)nSNPCnt/(float)nTotalRD;
        fReadDenScore = (fRS > 0.25)?100:fRS*400;

        if (nMappableReptDens==0) 
           nReptDenScore =20;
        else 
           nReptDenScore = ((float)nMappableReptDens/(float)nMappableBases < 66.00)?20:0; //16+(0.5*100) => 66.00
        
        fAQSScore = (fAvgQS < 15.00)?0:(fAvgQS-15.00)*(100.00/15.00);
        fSNPScore = (fPEScore+fSNPReadScore+fStrandScore+fReadDenScore+(float)nReptDenScore+fAQSScore)/520.00;
        
        fprintf(g_pfOut,"%s\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.4f\n",
                        acOut,fPEScore,fSNPReadScore,fStrandScore,fReadDenScore,nReptDenScore,fAQSScore,fSNPScore);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv) {

    if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_fMean = atof(argv[2]);

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

