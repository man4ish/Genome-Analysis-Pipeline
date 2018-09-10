/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on Dec 27, 2010, 10:30 AM
 */

#include "SXMergePoly.h"

/*
 * 
 */

FILE *g_pfIn, *g_pfOut;
char g_cSNVType;


float CalcAvg(float fTotalQS, float fTotalCnt)
{
    if (fTotalCnt == 0.00) return 0;
    return (fTotalQS == 0.00)? 0.00:fTotalQS/fTotalCnt;
}

void GenerateSNPList()
{
    fprintf(stdout,"GenerateSNPList\n"); 
    char acbuf[4096],acOut[4096]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    SNP_LIST::iterator itr; int nChromosome = 0; unsigned long ulTotalCnt;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateSNPList-Total Recs Processed = %lu M\n", ulTest2);
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

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);
        
        pstData = new stData;
        pChr = strtok(NULL,"\t"); pstData->nFwdCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nRvsCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nReadCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->unReadDen = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nPEndCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->fAvgQryPos = atof(pChr);
        pChr = strtok(NULL,"\t"); pstData->fAvgQScore = atof(pChr);
        pChr = strtok(NULL,"\t"); pstData->pcdbSNP = new char[strlen(pChr)+1];
        strcpy(pstData->pcdbSNP,pChr);
        
        itr = g_aSNPList[nChromosome-1].find(pstKey);  

        if (itr == g_aSNPList[nChromosome-1].end()){
           g_aSNPList[nChromosome-1][pstKey]=pstData; //fprintf(g_pfOut,"%s\n",acOut); continue;             
        } 
        else
        {
             ulTotalCnt = (*itr->second).nReadCnt + pstData->nReadCnt;               

             (*itr->second).fAvgQryPos = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQryPos)+
                                                 (pstData->nReadCnt * pstData->fAvgQryPos),ulTotalCnt); 
                                          
             (*itr->second).fAvgQScore = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQScore)+
                                                 (pstData->nReadCnt * pstData->fAvgQScore),ulTotalCnt);

             (*itr->second).nFwdCnt += pstData->nFwdCnt;
             (*itr->second).nRvsCnt += pstData->nRvsCnt;
             (*itr->second).nReadCnt += pstData->nReadCnt;
             (*itr->second).nPEndCnt += pstData->nPEndCnt;

             if (g_cSNVType == 'D'){
                if ((*itr->second).unReadDen > (unsigned int)pstData->nReadCnt) 
                   (*itr->second).unReadDen -= (unsigned int)pstData->nReadCnt; 
                else 
                   (*itr->second).unReadDen =0;
             }                     
        }  
    }
}


void OutputData()
{
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    SNP_LIST::iterator itr;
 
    for (int i=0; i <24; i++){
        itr = g_aSNPList[i].begin();

        while (itr != g_aSNPList[i].end())
        {
            ulTest1++;
            if (ulTest1 == 1000000)
            {
               ulTest2++;
               fprintf(stderr,"Output-Total Recs Processed = %lu M\n", ulTest2);
               ulTest1=0;
            }

            if (i < 22) fprintf(g_pfOut,"%d\t",i+1);
            else if (i == 22) fprintf(g_pfOut,"X\t");
            else fprintf(g_pfOut,"Y\t");

            fprintf(g_pfOut,"%u\t%c\t%u\t%u\t%u\t%u\t%u\t%.2f\t%.2f\t%s\n",
                             (*itr->first).nOffset,(*itr->first).acAllele[0],(*itr->second).nFwdCnt,
                             (*itr->second).nRvsCnt,(*itr->second).nReadCnt,(*itr->second).unReadDen,
                             (*itr->second).nPEndCnt,(*itr->second).fAvgQryPos,(*itr->second).fAvgQScore,
                             (*itr->second).pcdbSNP);
            itr++;
        }
    }    
}


int main(int argc, char** argv) {

    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");  
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input filename> <SNVType: I-INS, D-DEL> <Output filename>\n\n");
       exit(9);   
    }   

    g_pfIn = fopen(argv[1],"r"); if (!g_pfIn) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}
    g_pfOut = fopen(argv[3],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...",argv[3]); goto ExitRtn;}
    g_cSNVType = argv[2][0];

    GenerateSNPList(); OutputData();

 ExitRtn:
    if (g_pfIn) fclose (g_pfIn); 
    if (g_pfOut) fclose (g_pfOut);

    for (int i=0; i<24;i++) 
         g_aSNPList[i].clear();

    return (EXIT_SUCCESS);
}

