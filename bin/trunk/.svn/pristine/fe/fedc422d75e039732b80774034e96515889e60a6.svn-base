#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL, *g_pfOut=NULL;


static void banner(char *argv[])
{   
    printf("Synamatix SXAppendScoreNS Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Output filename>\n\n");
}


void ProcessRecs()
{
    int nFwdSNPCnt,nRvsSNPCnt,nSNPCnt,nPECnt,nTotalRD,nMappableBases,nMappableReptDens;
    int nPEScore,nSNPReadScore,nStrandScore,nReptDenScore;  
    float fAvgQS,fRS,fReadDenScore,fAQSScore,fSNPScore;  

    char acbuf[40960],acOut[40960],*pChr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {        
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
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

        pChr = strtok(acbuf,"\t"); //Chromosome
        pChr = strtok(NULL,"\t");  pChr = strtok(NULL,"\t"); //Offset,Nucleotide_Variant

        pChr = strtok(NULL,"\t");  nFwdSNPCnt = atoi(pChr); //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  nRvsSNPCnt = atoi(pChr); //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  nSNPCnt = atoi(pChr);    //SNP_Reads                 
        pChr = strtok(NULL,"\t");  nTotalRD = atoi(pChr);   //Total_Read_Density
        pChr = strtok(NULL,"\t");  nPECnt = atoi(pChr);     //PEnd_Count 

        pChr = strtok(NULL,"\t"); //Avg_QryPos 
        pChr = strtok(NULL,"\t"); fAvgQS = atof(pChr);//Avg_QScore
        //pChr = strtok(NULL,"\t"); //dbSNP      

        //pChr = strtok(NULL,"\t"); //Local_Copy_Number
        
        //pChr = strtok(NULL,"\t"); nMappableBases = atoi(pChr);//Mappable_Bases      
        //pChr = strtok(NULL,"\t"); nMappableReptDens = atoi(pChr);//Mappable_Repeat_Densities

        nPEScore = (nPECnt > 4)?100:nPECnt*20;
        nSNPReadScore = (nSNPCnt > 9)?100:nSNPCnt*10;
        nStrandScore = ((nFwdSNPCnt > 1) && (nRvsSNPCnt > 1))?100:0;  
        fRS = (float)nSNPCnt/(float)nTotalRD;
        fReadDenScore = (fRS > 0.25)?100:fRS*400;

        //if (nMappableReptDens==0) 
           nReptDenScore =20; //cause nMappableReptDens is always 0 for mouse sample
        //else 
        //   nReptDenScore = ((float)nMappableReptDens/(float)nMappableBases < 66.00)?20:0; //16+(0.5*100) => 66.00
        
        fAQSScore = (fAvgQS < 15.00)?0:(fAvgQS-15.00)*100/25;
        fSNPScore = ((float)nPEScore+(float)nSNPReadScore+(float)nStrandScore+fReadDenScore+(float)nReptDenScore+fAQSScore)/400.00;
        
        fprintf(g_pfOut,"%s\t%d\t%d\t%d\t%.2f\t%d\t%.2f\t%.4f\n",
                        acOut,nPEScore,nSNPReadScore,nStrandScore,fReadDenScore,nReptDenScore,fAQSScore,fSNPScore);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv) {

    if (argc < 3){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[2],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

