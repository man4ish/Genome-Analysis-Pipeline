#include "SXMergeSNPSample2INDEL.h"

FILE *g_pfInSNP, *g_pfInGapList, *g_pfInINDEL, *g_pfOut;
char g_cSNVType, acFOutName[1024]; int g_nSampleNo=0;


float CalcAvg(float fTotalQS, float fTotalCnt)
{
    if (fTotalCnt == 0.00) return 0;
    return (fTotalQS == 0.00)? 0.00:fTotalQS/fTotalCnt;
}


void GenerateSNPList()
{
    fprintf(stdout,"GenerateSNPList\n");
    char acbuf[4096],acOut[4096]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    SNP_LIST::iterator itr; int nChromosome=0; unsigned long ulTotalCnt;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfInSNP))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateSNPList-Total Recs Processed = %lu M\n", ulTest2);
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

        /*
        pChr = strtok(NULL,"\t");  if (g_nSampleNo == 1) {pstData->nReadCnt = atoi(pChr);}
        pChr = strtok(NULL,"\t");  if (g_nSampleNo == 2) {pstData->nReadCnt = atoi(pChr);} 
        pChr = strtok(NULL,"\t"); //pstData->nReadCnt = atoi(pChr); //Total_SNP_Reads

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstData->unReadDen = atoi(pChr); //Total_Read_Density

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstData->fAvgQScore=atof(pChr);//Total_Avg_QScore

        pChr = strtok(NULL,"\t"); //dbSnp
        pChr = strtok(NULL,"\t"); // Read Strength

        pChr = strtok(NULL,"\t"); if (g_nSampleNo == 1) {pstData->nFwdCnt = atoi(pChr);}
        pChr = strtok(NULL,"\t"); if (g_nSampleNo == 2) {pstData->nFwdCnt = atoi(pChr);}
        pChr = strtok(NULL,"\t"); //pstData->nFwdCnt = atoi(pChr); //Total_Fwd_SNP

        pChr = strtok(NULL,"\t"); if (g_nSampleNo == 1) {pstData->nRvsCnt = atoi(pChr);}
        pChr = strtok(NULL,"\t"); if (g_nSampleNo == 2) {pstData->nRvsCnt = atoi(pChr);}
        pChr = strtok(NULL,"\t"); //pstData->nRvsCnt = atoi(pChr); //Total_Rvs_SNP

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstData->nPEndCnt = atoi(pChr); //Total_PEnd

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstData->fAvgQryPos= atof(pChr); //Total_AVQryPos
        */ 

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
                if ((*itr->second).unReadDen > (unsigned int)pstData->nReadCnt)
                   (*itr->second).unReadDen -= (unsigned int)pstData->nReadCnt;
                else
                   (*itr->second).unReadDen =0;
             }
        }
    }
/*
    for (int i=0; i<24; i++)
    {
        fprintf(stdout,"Chro-%d\t%d\n", i+1, g_aSNPList[i].size());
    }
*/
}


void GenerateGapList()
{
    fprintf(stdout,"GenerateGapList\n");
    char acbuf[4096],acInfo[1024],acSNPInfo[1024],*pChr=NULL;
    stKey *pstKey=NULL; stGapData *pstGapData=NULL;
    int nChromosome=0;
    GAP_LIST::iterator itr;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfInGapList))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stderr,"GenerateGapList-Total Recs Processed = %lu M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInGapList)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;        
        
        pstKey = new stKey; pstGapData = new stGapData;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); strcpy(acSNPInfo,pChr);
        pChr = strtok(NULL,"\t");

        if (g_cSNVType == 'I') strcpy(acInfo,pChr);
        else if (g_cSNVType == 'D'){ pChr = strtok(NULL,"\t"); strcpy(acInfo,pChr); }

        pChr = strtok(NULL,"\t"); pstGapData->nGap=atoi(pChr);

        pChr = strtok(acInfo,":"); pstKey->nOffset = atoi(pChr+2);

        //itr = g_aGapList[nChromosome-1].find(pstKey);
        //if (itr == g_aGapList[nChromosome-1].end())

        pChr = strtok(acSNPInfo,":"); //nSNPOfs = atoi(pChr+2); 

        pstGapData->nSNPOfs = atoi(pChr+2);
        g_aGapList[nChromosome-1][pstKey]= pstGapData;//nSNPOfs;
    }
}


void OutputData()
{
    char acbuf[4096],acOut[4096], acChro[3], cAllele, *pChr=NULL;
    int nChromosome = 0, nOffset; unsigned long ulTotalCnt;
    stKey *pstKey=new stKey; stData *pstData = new stData; SNP_LIST::iterator itr;
    GAP_LIST::iterator gap_itr;

   int nCnt1=0,nCnt2=0, ndbSnp=0, nNovel=0; 
/*
char acTestName[1024];

sprintf(acTestName,"%s_Novel_SR_Added_%c",acFOutName,g_cSNVType);
FILE *pfNovel = fopen(acTestName,"w");

sprintf(acTestName,"%s_dbSnp_SR_Added_%c",acFOutName,g_cSNVType);
FILE *pfSnp = fopen(acTestName,"w");

FILE *pfTest=fopen("SNPOfs.log","w");
*/
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
    while(!feof(g_pfInINDEL))
    {
        if (++ulTest1 == 1000000)
        { fprintf(stderr,"OutputData-Total Recs Processed = %lu M\n", ++ulTest2); ulTest1=0; }

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

        pstKey->nOffset = (*gap_itr->second).nSNPOfs; //gap_itr->second; 
//fprintf(pfTest,"%d\t%u\n",nChromosome,pstKey->nOffset);

        itr = g_aSNPList[nChromosome-1].find(pstKey);

        if (itr == g_aSNPList[nChromosome-1].end()) { fprintf(g_pfOut,"%s\n",acOut); nCnt1++;}
        else{
           nCnt2++;
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

           ulTotalCnt = pstData->nReadCnt + (*itr->second).nReadCnt;
           pstData->fAvgQryPos = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQryPos)+
                                                 (pstData->nReadCnt * pstData->fAvgQryPos),ulTotalCnt);

           pstData->fAvgQScore = CalcAvg(((*itr->second).nReadCnt * (*itr->second).fAvgQScore)+
                                                (pstData->nReadCnt * pstData->fAvgQScore),ulTotalCnt);

           pstData->nFwdCnt += (*itr->second).nFwdCnt;
           pstData->nRvsCnt += (*itr->second).nRvsCnt;
           pstData->nPEndCnt += (*itr->second).nPEndCnt;

           if (g_cSNVType == 'D'){
                if (pstData->unReadDen > (unsigned int)((*itr->second).nReadCnt))
                    pstData->unReadDen -= (unsigned int)((*itr->second).nReadCnt);
                else
                    pstData->unReadDen = 0;
           }
   
           fprintf(g_pfOut,"%s\t%u\t%c\t%u\t%u\t%lu\t%u\t%u\t%.2f\t%.2f\t%s\n",
                            acChro,nOffset,cAllele,pstData->nFwdCnt,
                            pstData->nRvsCnt,ulTotalCnt,pstData->unReadDen,
                            pstData->nPEndCnt,pstData->fAvgQryPos,pstData->fAvgQScore,
                            pChr);//pstData->pcdbSNP);

//Testing ////////////////////////////////////////////////////////////////////////////////////////////
/*
if (pChr[0]=='-') {    
    nNovel++; 
    fprintf(pfNovel,"%s\t%u\t%c\t%u\t%u\t%u\t%u\t%u\t%.2f\t%.2f\t%s\t%u\t%u\t%u\t%u\n",
                            acChro,nOffset,cAllele,(*itr->second).nFwdCnt,
                            (*itr->second).nRvsCnt,ulTotalCnt,pstData->unReadDen,
                            (*itr->second).nPEndCnt,(*itr->second).fAvgQryPos,(*itr->second).fAvgQScore,
                            pChr,(*itr->second).nReadCnt,pstData->nReadCnt,pstKey->nOffset,(*gap_itr->second).nGap);
}
else {
    ndbSnp++;
    fprintf(pfSnp,"%s\t%u\t%c\t%u\t%u\t%u\t%u\t%u\t%.2f\t%.2f\t%s\t%u\t%u\t%u\t%u\n",
                            acChro,nOffset,cAllele,(*itr->second).nFwdCnt,
                            (*itr->second).nRvsCnt,ulTotalCnt,pstData->unReadDen,
                            (*itr->second).nPEndCnt,(*itr->second).fAvgQryPos,(*itr->second).fAvgQScore,
                            pChr,(*itr->second).nReadCnt,pstData->nReadCnt,pstKey->nOffset,(*gap_itr->second).nGap);
}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////

            itr++;
        }
    }

    if (pstKey) delete pstKey; if (pstData) delete pstData;

    fprintf(stdout,"nCnt1 => %d\n",nCnt1); fprintf(stdout,"nCnt2 => %d\n",nCnt2);  
    fprintf(stdout,"ndbSnp => %d\n",ndbSnp); fprintf(stdout,"nNovel => %d\n",nNovel);
}


int main(int argc, char *argv[])
{
    if (argc < 7)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input SNP filtered list> <Sample No.> <Input filename-records with gap specified> "
                      "<Input INS/DEL filename> <SNVType: I-INS, D-DEL> <Output filename>\n\n");
       exit(9);
    }

    g_pfInSNP = fopen(argv[1],"r"); if (!g_pfInSNP) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto ExitRtn;}
    g_nSampleNo = atoi(argv[2]); 
    g_pfInGapList = fopen(argv[3],"r"); if (!g_pfInGapList) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto ExitRtn;} 
    g_pfInINDEL = fopen(argv[4],"r"); if (!g_pfInINDEL) {fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto ExitRtn;}     
    g_cSNVType = argv[5][0];
    g_pfOut = fopen(argv[6],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[6]); goto ExitRtn;}
    strcpy(acFOutName,argv[6]);

    GenerateSNPList(); GenerateGapList(); OutputData();

 ExitRtn:
    if (g_pfInSNP) fclose (g_pfInSNP); if (g_pfInINDEL) fclose (g_pfInINDEL);
    if (g_pfOut) fclose (g_pfOut);

    for (int i=0; i<24;i++)
         g_aSNPList[i].clear();

    return (EXIT_SUCCESS);

}
