/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on May 25, 2010, 10:30 AM
 */

#include "SXMergePoly.h"

/*
 * 
 */

SNP_LIST g_SNPList;
FILE *g_pfIn, *g_pfOut;


void GenerateSNPList()
{
    fprintf(stdout,"GenerateSNPList\n");
    char acbuf[4096],acOut[4096]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    SNP_LIST::iterator itr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateSNPList-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
        {
           fprintf(g_pfOut,"%s\n",acbuf); continue;
        }
 
        strcpy(acOut,acbuf); pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) pstKey->cChromosome = atoi(pChr);
        else if (pChr[0]=='X') pstKey->cChromosome = 23;
        else if(pChr[0]=='Y') pstKey->cChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);
        
        pstData = new stData;
        pChr = strtok(NULL,"\t"); pstData->nSNPCnt = atoi(pChr);
        pChr = strtok(NULL,"\t");  pstData->nReadDensity = atoi(pChr);
        pChr = strtok(NULL,"\t");  pstData->fAvgQScore = atof(pChr);
        pChr = strtok(NULL,"\t"); pstData->pcdbSNP = new char[strlen(pChr)+1];
        strcpy(pstData->pcdbSNP,pChr);
        
        itr = g_SNPList.find(pstKey);  
        if (itr == g_SNPList.end()){
           g_SNPList[pstKey]=pstData; fprintf(g_pfOut,"%s\n",acOut); continue;  
           
        } 
    }
}


int main(int argc, char** argv) {

    if (argc < 3)
    {
       fprintf(stdout,"Invalid Parameters...\n");  
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input filename> <Output filename> ");
       exit(9);   
    }   

    g_pfIn = fopen(argv[1],"r"); if (!g_pfIn) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}
    g_pfOut = fopen(argv[2],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...",argv[2]); goto ExitRtn;}

    GenerateSNPList();

 ExitRtn:
    if (g_pfIn) fclose (g_pfIn); 
    if (g_pfOut) fclose (g_pfOut);

    g_SNPList.clear();

    return (EXIT_SUCCESS);
}

