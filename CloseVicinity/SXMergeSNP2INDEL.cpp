#include "SXMergeSNP2INDEL.h"

FILE *g_pfInSNP, *g_pfInGapList, *g_pfInINDEL, *g_pfOut;
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
    SNP_LIST::iterator itr=NULL; int nChromosome; unsigned long ulTotalCnt;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfInSNP))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateSNPList-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInSNP)){break;}

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

        //pChr = strtok(NULL,"\t"); //giNumber
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); //strcpy(pstKey->acAllele,pChr);

        pstData = new stData;
        pChr = strtok(NULL,"\t"); pstData->nFwdCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nRvsCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nReadCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->unReadDen = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nPEndCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->fAvgQryPos = atof(pChr);
        pChr = strtok(NULL,"\t"); pstData->fAvgQScore = atof(pChr);
        //pChr = strtok(NULL,"\t"); pstData->pcdbSNP = new char[strlen(pChr)+1];
        //strcpy(pstData->pcdbSNP,pChr);

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
                if ((*itr->second).unReadDen > pstData->nReadCnt)
                   (*itr->second).unReadDen -= pstData->nReadCnt;
                else
                   (*itr->second).unReadDen =0;
             }
        }
    }
/*
    for (int i=0; i<24; i++)
    {
        fprintf(stdout,"Chro-%d\n", i+1, g_aSNPList[i].size());
    }
*/
}


void GenerateGapList()
{
    fprintf(stdout,"GenerateGapList\n");
    char acbuf[4096],acInfo[1024],acSNPInfo[1024],*pChr=NULL,cAllele;
    stKey *pstKey=NULL;
    int nChromosome, nSReads, nSNPOfs; unsigned long ulTotalCnt;
    GAP_LIST::iterator itr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfInGapList))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stderr,"GenerateGapList-Total Recs Processed = %u M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInGapList)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;        
        
        pstKey = new stKey; pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); strcpy(acSNPInfo,pChr);
        pChr = strtok(NULL,"\t");

        if (g_cSNVType == 'I') strcpy(acInfo,pChr);
        else if (g_cSNVType == 'D'){ pChr = strtok(NULL,"\t"); strcpy(acInfo,pChr); }

        pChr = strtok(acInfo,":"); pstKey->nOffset = atoi(pChr+2);

        //itr = g_aGapList[nChromosome-1].find(pstKey);
        //if (itr == g_aGapList[nChromosome-1].end())

        pChr = strtok(acSNPInfo,":"); nSNPOfs = atoi(pChr+2); 

        g_aGapList[nChromosome-1][pstKey]=nSNPOfs;
    }
}


void OutputData()
{
    char acbuf[4096],acOut[4096],acInfo[1024], acChro[3], cAllele, *pChr=NULL;
    int nChromosome, nOffset, nSR, nDiff; unsigned long ulTotalCnt;
    stKey *pstKey=new stKey; stData *pstData = new stData; SNP_LIST::iterator itr=NULL;
    GAP_LIST::iterator gap_itr=NULL;
    
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
    while(!feof(g_pfInINDEL))
    {
        if (++ulTest1 == 1000000)
        { fprintf(stderr,"OutputData-Total Recs Processed = %u M\n", ++ulTest2); ulTest1=0; }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInINDEL)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0]==0||acbuf[0]=='C'||acbuf[0]=='#') continue;

        strcpy(acOut,acbuf);         

        pChr = strtok(acbuf,"\t"); strcpy(acChro,pChr);

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t");//gi_Number
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr); pstKey->nOffset = nOffset;

        gap_itr = g_aGapList[nChromosome-1].find(pstKey); 

        if (gap_itr == g_aGapList[nChromosome-1].end()) {fprintf(g_pfOut,"%s\n",acOut); continue;}

        pstKey->nOffset = gap_itr->second; 

        itr = g_aSNPList[nChromosome-1].find(pstKey);

        if (itr == g_aSNPList[nChromosome-1].end()) { fprintf(g_pfOut,"%s\n",acOut); }
        else{
           pChr = strtok(NULL,"\t"); cAllele = pChr[0];
           pChr = strtok(NULL,"\t"); pstData->nFwdCnt = atoi(pChr);
           pChr = strtok(NULL,"\t"); pstData->nRvsCnt = atoi(pChr);           
           pChr = strtok(NULL,"\t"); pstData->nReadCnt = atoi(pChr);

           //if (pstData->nReadCnt < (*itr->second).nReadCnt) {fprintf(g_pfOut,"%s\n",acOut); continue;}

           pChr = strtok(NULL,"\t"); pstData->unReadDen = atoi(pChr);
           pChr = strtok(NULL,"\t"); pstData->nPEndCnt = atoi(pChr);
           pChr = strtok(NULL,"\t"); pstData->fAvgQryPos = atof(pChr);
           pChr = strtok(NULL,"\t"); pstData->fAvgQScore = atof(pChr);
           pChr = strtok(NULL,"\t"); //dbSnp

           //if (pChr[0]=='-') {fprintf(g_pfOut,"%s\n",acOut); continue;} //only merge the known dbSnp rec.

           //pstData->pcdbSNP = new char[strlen(pChr)+1];
           //strcpy(pstData->pcdbSNP,pChr);  

           ulTotalCnt = (*itr->second).nReadCnt + pstData->nReadCnt;
           (*itr->second).fAvgQryPos = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQryPos)+
                                                 (pstData->nReadCnt * pstData->fAvgQryPos),ulTotalCnt);

           (*itr->second).fAvgQScore = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQScore)+
                                                (pstData->nReadCnt * pstData->fAvgQScore),ulTotalCnt);

           (*itr->second).nFwdCnt += pstData->nFwdCnt;
           (*itr->second).nRvsCnt += pstData->nRvsCnt;
           (*itr->second).nPEndCnt += pstData->nPEndCnt;

           if (g_cSNVType == 'D'){
                if (pstData->unReadDen > (*itr->second).nReadCnt)
                    pstData->unReadDen -= (*itr->second).nReadCnt;
                else
                    pstData->unReadDen = 0;
  
                /*if ((*itr->second).unReadDen > pstData->nReadCnt)
                   (*itr->second).unReadDen -= pstData->nReadCnt;
                else
                   (*itr->second).unReadDen =0;
                */
           }
   
           fprintf(g_pfOut,"%s\t%u\t%c\t%u\t%u\t%u\t%u\t%u\t%.2f\t%.2f\t%s\n",
                            acChro,nOffset,cAllele,(*itr->second).nFwdCnt,
                            (*itr->second).nRvsCnt,ulTotalCnt,pstData->unReadDen,
                            (*itr->second).nPEndCnt,(*itr->second).fAvgQryPos,(*itr->second).fAvgQScore,
                            pChr);//pstData->pcdbSNP);
            itr++;
        }
    }

    if (pstKey) delete pstKey; if (pstData) delete pstData;
}


int main(int argc, char *argv[])
{
    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input SNP filtered list> <Input filename-records with gap specified> "
                      "<Input INS/DEL filename> <SNVType: I-INS, D-DEL> <Output filename>\n\n");
       exit(9);
    }

    g_pfInSNP = fopen(argv[1],"r"); if (!g_pfInSNP) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto ExitRtn;}
    g_pfInGapList = fopen(argv[2],"r"); if (!g_pfInGapList) {fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto ExitRtn;} 
    g_pfInINDEL = fopen(argv[3],"r"); if (!g_pfInINDEL) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto ExitRtn;}     
    g_cSNVType = argv[4][0];
    g_pfOut = fopen(argv[5],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[5]); goto ExitRtn;}

    GenerateSNPList(); GenerateGapList(); OutputData();

 ExitRtn:
    if (g_pfInSNP) fclose (g_pfInSNP); if (g_pfInINDEL) fclose (g_pfInINDEL);
    if (g_pfOut) fclose (g_pfOut);

    for (int i=0; i<24;i++)
         g_aSNPList[i].clear();

    return (EXIT_SUCCESS);

}
