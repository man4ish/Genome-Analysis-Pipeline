#include "SXMergeV2.h"

stRecList g_pRecs[24];
FILE *g_pfFile1, *g_pfFile2, *g_pfCommon, *g_pfUnique1, *g_pfUnique2;


void GenerateList()
{
    fprintf(stdout,"Generating Data List...\n");
    char acbuf[1024]; char *pChr=NULL; int nChromosome; stRecList *pRecList=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfFile1))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Generating Data List-Total Recs = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFile1)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
           continue;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pRecList = &g_pRecs[nChromosome-1];   

        if (pRecList->nUsed == pRecList->nTotal){
            pRecList->nTotal+=LIST_SIZE;
            pRecList->pData = (stData*)realloc(pRecList->pData,sizeof(stData)*pRecList->nTotal);
        }

        pRecList->pData[pRecList->nUsed].bUnique=true;  
        pChr = strtok(NULL,"\t"); pRecList->pData[pRecList->nUsed].nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pRecList->pData[pRecList->nUsed].acAllele,pChr);
        pChr = strtok(NULL,"\t"); pRecList->pData[pRecList->nUsed].nSNPCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); pRecList->pData[pRecList->nUsed].nReadDensity = atoi(pChr);
        pChr = strtok(NULL,"\t"); pRecList->pData[pRecList->nUsed].fAvgQScore = atof(pChr);
        pChr = strtok(NULL,"\t"); pRecList->pData[pRecList->nUsed].pcdbSNP = new char[strlen(pChr)+1];
        strcpy(pRecList->pData[pRecList->nUsed].pcdbSNP,pChr); pRecList->nUsed++;        
    }
}


int SearchData(stRecList *pRecList, stData &ostData)
{
    int nLow=0,nHigh=pRecList->nUsed,nPos,nIdx=-1;       
    
    while (!(nLow > nHigh))
    {
    	nPos = (nHigh+nLow)/2;
        if (pRecList->pData[nPos].nOffset == ostData.nOffset)
        {
            for (int i=nPos-2; i < nPos+3; i++)
            {
                 if (i < 0) continue; 
                 if (pRecList->pData[i].nOffset != ostData.nOffset) continue;
                 if (strcmp(pRecList->pData[i].acAllele,ostData.acAllele) != 0) continue;
                 nIdx = i; break;
            }   
            break;           
        } 
        else if (pRecList->pData[nPos].nOffset < ostData.nOffset) nLow = nPos+1;
        else nHigh = nPos-1;                  
    } 

    return nIdx;
}


void AnalyseData()
{
    fprintf(stdout,"Analysing Data ... \n");
    char acbuf[1024]; char *pChr=NULL; stRecList *pRecList=NULL; stData ostData;
    unsigned long ulTest1=0,ulTest2=0; int nChromosome,nReadStrength,nIdx;         

    while(!feof(g_pfFile2))
    {
        ulTest1++;
        if (ulTest1 == 1000000){
            ulTest2++;
            fprintf(stderr,"Analysing Data-Total Recs = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFile2)){break;}      
    
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pRecList=&g_pRecs[nChromosome-1];

        pChr = strtok(NULL,"\t"); ostData.nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(ostData.acAllele,pChr);
        pChr = strtok(NULL,"\t"); ostData.nSNPCnt = atoi(pChr);
        pChr = strtok(NULL,"\t"); ostData.nReadDensity = atoi(pChr);
        pChr = strtok(NULL,"\t"); ostData.fAvgQScore = atof(pChr);

        nIdx = SearchData(pRecList,ostData);  

        if (nIdx==-1)
        {              
            pChr = strtok(NULL,"\t"); //dbSNP

            if (ostData.nReadDensity < 0) nReadStrength = 101;
            else nReadStrength = (int)(ostData.nSNPCnt*100/ostData.nReadDensity+0.5);

            if (nChromosome < 23) fprintf(g_pfUnique2,"%d\t",nChromosome); 
            else if(nChromosome == 23) fprintf(g_pfUnique2,"X\t"); else fprintf(g_pfUnique2,"Y\t");

            fprintf(g_pfUnique2,"%d\t%s\t-\t%d\t%d\t-\t%d\t%d\t-\t%.2f\t%.2f\t%s\t%d\n",
                    ostData.nOffset,ostData.acAllele,ostData.nSNPCnt,ostData.nSNPCnt,
                    ostData.nReadDensity,ostData.nReadDensity,ostData.fAvgQScore,
                    ostData.fAvgQScore,pChr,nReadStrength);//(int)(pstData->nSNPCnt*100/pstData->nReadDensity+0.5));
        }
        else{
            pRecList->pData[nIdx].bUnique = false;
            pRecList->pData[nIdx].nSNPCnt2 = ostData.nSNPCnt;
            pRecList->pData[nIdx].nReadDensity2 = ostData.nReadDensity;
            pRecList->pData[nIdx].fAvgQScore2 = ostData.fAvgQScore;
        }
    }
}


void OutputData()
{
    fprintf(g_pfCommon,"Chr\tOffset\tAllele\tSNP_Density_1\tSNP_Density_2\t"
                       "Total_SNP_ Density\tRead_Density_1\tRead_Density_2\t"
                       "Total_Read_Density\tAVQS_1\tAVQS_2\tAVQS\tdbSNP\tRS\n");

    fprintf(g_pfUnique1,"Chr\tOffset\tAllele\tSNP_Density_1\tSNP_Density_2\t"
                        "Total_SNP_ Density\tRead_Density_1\tRead_Density_2\t"
                        "Total_Read_Density\tAVQS_1\tAVQS_2\tAVQS\tdbSNP\tRS\n");

    fprintf(g_pfUnique2,"Chr\tOffset\tAllele\tSNP_Density_1\tSNP_Density_2\t"
                        "Total_SNP_ Density\tRead_Density_1\tRead_Density_2\t"
                        "Total_Read_Density\tAVQS_1\tAVQS_2\tAVQS\tdbSNP\tRS\n");

    fprintf(stdout,"Output Data...\n");
    int nTotalSNPCnt, nTotalReadDensity, nReadStrength;
    FILE *pf=NULL; stRecList *pstRecList=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    for (int i=0; i < 24; i++)
    {
       pstRecList = &g_pRecs[i];

       for (int j=0; j< pstRecList->nUsed; j++)        
       {
            ulTest1++;
            if (ulTest1 == 1000000){
                ulTest2++;
                fprintf(stderr,"Total Recs Processed = %u M\n", ulTest2);
                ulTest1=0;
            }

            pf = (pstRecList->pData[j].bUnique)? g_pfUnique1 : g_pfCommon;
  
            if (i < 22) fprintf(pf,"%d\t",i+1); else if(i == 22) fprintf(pf,"X\t"); else fprintf(pf,"Y\t");

            if (pstRecList->pData[j].bUnique)
            {
                if (pstRecList->pData[j].nReadDensity==0) nReadStrength = 101;
                else nReadStrength = (int)(((float)pstRecList->pData[j].nSNPCnt*100/(float)pstRecList->pData[j].nReadDensity)+0.5);

                fprintf(g_pfUnique1,"%d\t%s\t%d\t-\t%d\t%d\t-\t%d\t%.2f\t-\t%.2f\t%s\t%d\n",
                        pstRecList->pData[j].nOffset,pstRecList->pData[j].acAllele,pstRecList->pData[j].nSNPCnt,
                        pstRecList->pData[j].nSNPCnt,pstRecList->pData[j].nReadDensity,
                        pstRecList->pData[j].nReadDensity,pstRecList->pData[j].fAvgQScore,
                        pstRecList->pData[j].fAvgQScore,pstRecList->pData[j].pcdbSNP,nReadStrength);       
            }
            else
            {
                nTotalSNPCnt = pstRecList->pData[j].nSNPCnt + pstRecList->pData[j].nSNPCnt2;
                nTotalReadDensity = pstRecList->pData[j].nReadDensity+pstRecList->pData[j].nReadDensity2;

                if (nTotalReadDensity == 0) nReadStrength = 101;
                else nReadStrength = (int)(((float)nTotalSNPCnt*100/(float)nTotalReadDensity)+0.5);

                fprintf(g_pfCommon,"%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%s\t%d\n",
                        pstRecList->pData[j].nOffset,pstRecList->pData[j].acAllele,
                        pstRecList->pData[j].nSNPCnt,pstRecList->pData[j].nSNPCnt2,nTotalSNPCnt,
                        pstRecList->pData[j].nReadDensity,pstRecList->pData[j].nReadDensity2,nTotalReadDensity,
                        pstRecList->pData[j].fAvgQScore,pstRecList->pData[j].fAvgQScore2,
                        (pstRecList->pData[j].fAvgQScore+pstRecList->pData[j].fAvgQScore2)/2,
                        pstRecList->pData[j].pcdbSNP,nReadStrength); 
            }
        }
       
    }
}



int main(int argc, char** argv) {

    /*
    <Input file name 1> <Input file name 2> <Output file name 1 - common data >
    <Output file name 2 - Unique data set 1> <Output file name 3 - Unique data set 2>
   */

    if (argc < 6)
    {
       fprintf(stdout,"Invalid Parameters...\n");  
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input filename 1> <Input filename 2> <Output filename 1 - common data > ");
       fprintf(stdout,"<Output filename 2 - Unique data set 1> <Output filename 3 - Unique data set 2>\n\n");
       exit(9);   
    }

    g_pfFile1=g_pfFile1=g_pfCommon=g_pfUnique1=g_pfUnique2=NULL;

    g_pfFile1 = fopen(argv[1],"r");
    if (!g_pfFile1) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}

    g_pfFile2 = fopen(argv[2],"r");
    if (!g_pfFile2) {fprintf(stdout,"Failed to open %s ...",argv[2]); goto ExitRtn;}

    g_pfCommon = fopen(argv[3],"w");
    if (!g_pfCommon) {fprintf(stdout,"Failed to open %s ...",argv[3]); goto ExitRtn;}
    g_pfUnique1 = fopen(argv[4],"w");
    if (!g_pfUnique1) {fprintf(stdout,"Failed to open %s ...",argv[4]); goto ExitRtn;}

    g_pfUnique2 = fopen(argv[5],"w");
    if (!g_pfUnique2) {fprintf(stdout,"Failed to open %s ...",argv[5]); goto ExitRtn;}

    GenerateList(); AnalyseData(); OutputData();

 ExitRtn:
    if (g_pfFile1) fclose (g_pfFile1);
    if (g_pfFile2) fclose (g_pfFile2);
    if (g_pfCommon) fclose (g_pfCommon);
    if (g_pfUnique1) fclose (g_pfUnique1);
    if (g_pfUnique2) fclose (g_pfUnique2);

    return (EXIT_SUCCESS);
}


