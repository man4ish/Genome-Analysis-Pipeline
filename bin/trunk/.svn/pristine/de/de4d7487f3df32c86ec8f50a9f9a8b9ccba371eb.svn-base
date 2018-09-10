#include <stdio.h>
#include <stdlib.h>
#include "GetCloseVicinitySNV.h"


void GenerateList(FILE *pfIn, SNV_LIST aSNVList[])
{
    fprintf(stdout,"GenerateList\n");
    char acbuf[4096], acOut[4096],*pChr=NULL, cAllele; 
    stKey *pstKey=NULL; 
    int nChromosome, nSReads; unsigned long ulTotalCnt;
    SNV_LIST::iterator itr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateList-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
        {
           continue;
        }

        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t"); //gi number
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
          
        itr = aSNVList[nChromosome-1].find(pstKey);

        if (itr==aSNVList[nChromosome-1].end())
        {
            aSNVList[nChromosome-1][pstKey]=1;
        } 
    }
}


void OutputData()
{
    fprintf(stdout,"OutputData\n");
    char acbuf[4096],acOut[4096], acChro[3], acAllele[4], *pChr=NULL; 
    stKey *pstKey=NULL; SNV_LIST::iterator itr=NULL; 
    int nChromosome, nOffset, nSR, nDiff; unsigned long ulTotalCnt;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; bool bINS_Found, bDEL_Found;

    while(!feof(g_pfIn_SNP))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"OutputData-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn_SNP)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
        {
           //fprintf(g_pfOut,"%s\n",acbuf); 
            continue;
        }

        strcpy(acOut,acbuf); pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        strcpy(acChro,pChr);
        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t"); // gi Number
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); // Nucleotide_Variant
        strcpy(acAllele,pChr);

        pChr = strtok(NULL,"\t"); // Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); // Rvs_SNP_Reads
        pChr = strtok(NULL,"\t"); // Total_SNP_Reads
        nSR = atoi(pChr); 

        bINS_Found = bDEL_Found = false;

        for (int i=nOffset-3; i < nOffset+4; i++)
        {
            pstKey->nOffset = i; 
            itr = g_aINSList[nChromosome-1].find(pstKey);

            if (itr!=g_aINSList[nChromosome-1].end())
            { 
                nDiff = abs(nOffset-i);                
                bINS_Found=true; break;
            }
        }

        if (bINS_Found) { fprintf(g_pfOut,"%s\n",acOut); continue; }
   
        for (int i=nOffset-3; i < nOffset+4; i++)
        {
            pstKey->nOffset = i; itr = g_aDELList[nChromosome-1].find(pstKey);

            if (itr!=g_aDELList[nChromosome-1].end())
            {
                nDiff = abs(nOffset-i);                                
                bDEL_Found=true; break;
            }
           
        }

        if (bDEL_Found) { fprintf(g_pfOut,"%s\n",acOut);continue; }
    }
}


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input filename-SNP List> <Input filename-INS List> <Input filename-DEL List> <Output filename>\n\n");
       exit(9);
    }
    
    g_pfIn_SNP = fopen(argv[1],"r"); if (!g_pfIn_SNP) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}
    g_pfIn_INS = fopen(argv[2],"r"); if (!g_pfIn_INS) {fprintf(stdout,"Failed to open %s ...",argv[2]); goto ExitRtn;}       
    g_pfIn_DEL = fopen(argv[3],"r"); if (!g_pfIn_DEL) {fprintf(stdout,"Failed to open %s ...",argv[3]); goto ExitRtn;}
    g_pfOut = fopen(argv[4],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...",argv[4]); goto ExitRtn;}

    GenerateList(g_pfIn_INS, g_aINSList); GenerateList(g_pfIn_DEL, g_aDELList); OutputData();


    //FILE *pfINSLog = fopen("INS.log","w");

    //SNV_LIST::iterator itr=NULL;
/*
    for (int i=0; i<24; i++)
    {
        itr = g_aINSList[i].begin();

        while (itr!=g_aINSList[i].end()) 
        {
            fprintf(pfINSLog,"%d\t%u\n",i+1,(*itr->first).nOffset);
            itr++;
        }
    }   
*/

 
 ExitRtn:
    if (g_pfIn_SNP) fclose (g_pfIn_SNP); if (g_pfIn_INS) fclose (g_pfIn_INS); 
    if (g_pfIn_DEL) fclose (g_pfIn_DEL); if (g_pfOut) fclose (g_pfOut);

    for (int i=0; i<24;i++) 
        g_aINSList[i].clear();

    return (EXIT_SUCCESS);

}
