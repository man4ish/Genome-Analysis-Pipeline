#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL, *g_pfOut=NULL; 
float g_fMean, g_fXCoverage, g_fYCoverage;


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
    int nStrandScore,nReptDenScore;  
    float fPEScore,fAvgQS,fRS,fReadDenScore,fAQSScore,fSNPScore,fSNPReadScore;  

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

        //fPEScore = !(nPECnt < g_fYCoverage)?100:(float)nPECnt*(100.00/g_fYCoverage); //(nPECnt > 2)?100:(float)nPECnt*(100.00/3.00);
        fSNPReadScore = !(nSNPCnt < g_fYCoverage)?100:nSNPCnt*(100.00/g_fYCoverage); 
        nStrandScore = (!(nFwdSNPCnt < g_fXCoverage) && !(nRvsSNPCnt < g_fXCoverage))?100:0;
        fRS = (float)nSNPCnt/(float)nTotalRD*100.00;
        fReadDenScore = (fRS > 25.00)?100.00:fRS*4.00; //(fRS > 10.00)?100.00:fRS*10.00;

        if (nMappableReptDens==0) 
           nReptDenScore =20;
        else 
           nReptDenScore = ((float)nMappableReptDens/(float)nMappableBases < 66.00)?20:0; //16+(0.5*100) => 66.00
        
        fAQSScore = (fAvgQS < 15.00)?0.00:(fAvgQS-15.00)*(100.00/25.00); //Denomitor = 40 (Max) - 15
        //fSNPScore = (fPEScore+fSNPReadScore+(float)nStrandScore+fReadDenScore+(float)nReptDenScore+fAQSScore)/520.00;

        fPEScore = 0.00;
        fSNPScore = (fSNPReadScore+(float)nStrandScore+fReadDenScore+(float)nReptDenScore+fAQSScore)/420.00;
        
        fprintf(g_pfOut,"%s\t%.2f\t%.2f\t%d\t%.2f\t%d\t%.2f\t%.4f\n",
                        acOut,fPEScore,fSNPReadScore,nStrandScore,fReadDenScore,nReptDenScore,fAQSScore,fSNPScore);
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

    g_fMean = atof(argv[2]); g_fXCoverage=g_fMean/5.00; g_fYCoverage=g_fMean/2.00;    

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

