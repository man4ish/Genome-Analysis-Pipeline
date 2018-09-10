#include <stdio.h>
#include <stdlib.h>
#include "SXChkSNVPoly.h"


void GenerateList(FILE *pfIn, SNV_LIST aSNVList[])
{
    fprintf(stdout,"GenerateList\n");
    char acbuf[40960],*pChr=NULL, cAllele; 
    stKey *pstKey=NULL; stData *pstData=NULL;
    int nChromosome, nSReads, nCnt; unsigned long ulTotalCnt;
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
           continue;

        //strcpy(acOut,acbuf); 
        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t"); //Gi Number
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
          
        pChr = strtok(NULL,"\t"); //Neucleotide_Variant
        cAllele = pChr[0];

        pChr = strtok(NULL,"\t"); //Fwd_INS_Reads
        pChr = strtok(NULL,"\t"); //Rvs_INS_Reads 
        pChr = strtok(NULL,"\t"); //Total_INS_Reads
         
        nSReads = atoi(pChr);

        itr = aSNVList[nChromosome-1].find(pstKey);

        if (itr==aSNVList[nChromosome-1].end())
        {
            pstData = new stData; 
            pstData->acAllele[0]=cAllele;
            pstData->anSReads[0]=nSReads;                    
            pstData->nCnt=1;
            aSNVList[nChromosome-1][pstKey]=pstData;
        } 
        else
        {
            (*itr->second).acAllele[(*itr->second).nCnt]=cAllele;
            (*itr->second).anSReads[1]=nSReads;      
            (*itr->second).nCnt++;
            
            delete pstKey;
        }
    }
}


void OutputData()
{
    fprintf(stdout,"OutputData\n");
    char acbuf[4096],acOut[4096], acChro[3], acAllele[4], *pChr=NULL; 
    stKey *pstKey=new stKey; SNV_LIST::iterator itr=NULL; 
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
            continue;

        strcpy(acOut,acbuf); 

        pChr = strtok(acbuf,"\t");

        strcpy(acChro,pChr);
        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t"); // Gi Number  
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); // Nucleotide_Variant
        strcpy(acAllele,pChr);

        pChr = strtok(NULL,"\t"); // Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); // Rvs_SNP_Reads
        pChr = strtok(NULL,"\t"); // Total_SNP_Reads
        nSR = atoi(pChr); 

        bINS_Found = bDEL_Found = false;

        fprintf(g_pfOut,"%s\tS-%u:%s,%u\tI-",acChro,nOffset,acAllele,nSR);
 
        for (int i=nOffset-3; i < nOffset+4; i++)
        {
            pstKey->nOffset = i; 
            itr = g_aINSList[nChromosome-1].find(pstKey);

            if (itr!=g_aINSList[nChromosome-1].end())
            { 
                nDiff = abs(nOffset-i);                
                fprintf(g_pfOut,"%u:",pstKey->nOffset);
                for (int j=0; j < (*itr->second).nCnt; j++)
                {
                    if (j > 0) fprintf(g_pfOut,"&");  
                    fprintf(g_pfOut,"%c,%u",(*itr->second).acAllele[j],(*itr->second).anSReads[j]);
                }   

                fprintf(g_pfOut,"\t%u\t",nDiff); bINS_Found=true; break;
            }
        }

        if (!bINS_Found) {fprintf(g_pfOut,"0\t");}

   
        fprintf(g_pfOut,"D-"); 
        for (int i=nOffset-3; i < nOffset+4; i++)
        {
            pstKey->nOffset = i; itr = g_aDELList[nChromosome-1].find(pstKey);

            if (itr!=g_aDELList[nChromosome-1].end())
            {
                nDiff = abs(nOffset-i);                                
                fprintf(g_pfOut,"%u:",pstKey->nOffset);

                for (int j=0; j < (*itr->second).nCnt; j++)
                {
                    if (j > 0) fprintf(g_pfOut,"&");
                    fprintf(g_pfOut,"%c,%u",(*itr->second).acAllele[j],(*itr->second).anSReads[j]);
                } 

                fprintf(g_pfOut,"\t%u\t",nDiff); bDEL_Found=true; break;
            }
           
        }

        if (!bDEL_Found) {fprintf(g_pfOut,"0");}  

        fprintf(g_pfOut,"\n"); 
    }

    if (pstKey) delete pstKey;
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
    
//FILE *pfINSLog = fopen("INS.log","w");
//SNV_LIST::iterator itr=NULL;

    g_pfIn_SNP = fopen(argv[1],"r"); if (!g_pfIn_SNP) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}
    g_pfIn_INS = fopen(argv[2],"r"); if (!g_pfIn_INS) {fprintf(stdout,"Failed to open %s ...",argv[2]); goto ExitRtn;}       
    g_pfIn_DEL = fopen(argv[3],"r"); if (!g_pfIn_DEL) {fprintf(stdout,"Failed to open %s ...",argv[3]); goto ExitRtn;}
    g_pfOut = fopen(argv[4],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...",argv[4]); goto ExitRtn;}

    GenerateList(g_pfIn_INS, g_aINSList); GenerateList(g_pfIn_DEL, g_aDELList); 
/*
    for (int i=0; i<24; i++)
        fprintf(stdout,"INS - %u\n",g_aINSList[i].size());

    for (int i=0; i<24; i++)
        fprintf(stdout,"DEL - %u\n",g_aDELList[i].size());
*/ 
    OutputData();

    fprintf(stdout,"Done...\n"); 

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

    for (int i=0; i<24;i++){ 
//fprintf(stdout,"Clear INS - %d\n",i+1);
        g_aINSList[i].clear(); 
//fprintf(stdout,"Clear DEL - %d\n",i+1);
        g_aDELList[i].clear();
     } 

    return (EXIT_SUCCESS);

}
