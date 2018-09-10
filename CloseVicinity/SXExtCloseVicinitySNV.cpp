#include <stdio.h>
#include <stdlib.h>
#include "SXExtCloseVicinitySNV.h"


void GenerateList(FILE *pfIn, SNV_LIST aSNVList[])
{
    fprintf(stdout,"GenerateList\n");
    char acbuf[40960],*pChr=NULL; 
    stKey *pstKey=NULL; 
    int nChromosome = 0, nSReads;
    SNV_LIST::iterator itr;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateList-Total Recs Processed = %lu M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
        {
           continue;
        }

        //strcpy(acOut,acbuf); 
        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
          
        pChr = strtok(NULL,"\t"); //Neucleotide_Variant
        //cAllele = pChr[0]; 
        pstKey->cAllele = pChr[0];

        pChr = strtok(NULL,"\t"); //Fwd_INS_Reads
        pChr = strtok(NULL,"\t"); //Rvs_INS_Reads 
        pChr = strtok(NULL,"\t"); //Total_INS_Reads
         
        nSReads = atoi(pChr);

        //itr = aSNVList[nChromosome-1].find(pstKey);

        //if (itr==aSNVList[nChromosome-1].end())
        //{
            /*
            pstData = new stData;
            pstData->acAllele[0]=cAllele;
            pstData->anSReads[0]=nSReads;                    
            pstData->nCnt=1;
            aSNVList[nChromosome-1][pstKey]=pstData;
            */

            aSNVList[nChromosome-1][pstKey]=nSReads;
        //} 
        //else
        //{
        //    (*itr->second).acAllele[(*itr->second).nCnt]=cAllele;
         //   (*itr->second).anSReads[(*itr->second).nCnt]=nSReads;      
         //   (*itr->second).nCnt++;

         //   delete pstKey;
        //}
    }
}


void OutputData()
{
    fprintf(stdout,"OutputData\n");
    char acbuf[4096],acOut[4096], acChro[3], acAllele[4], *pChr=NULL; 
    stKey *pstKey=new stKey; 
    SNV_LIST::iterator itr_INS; SNV_LIST::iterator itr_DEL;
    int nChromosome = 0, nOffset, nOffset_INS=0, nOffset_DEL=0, nSR, nDiff_INS=0, nDiff_DEL=0;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; bool bINS_Found, bDEL_Found;

    while(!feof(g_pfIn_SNP))
    {
        
        if (++ulTest1==1000000){fprintf(stderr,"OutputData-Total Recs Processed = %lu M\n", ++ulTest2);ulTest1=0;}

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

        //pChr = strtok(NULL,"\t"); // gi Number
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); // Nucleotide_Variant
        strcpy(acAllele,pChr);

        pChr = strtok(NULL,"\t"); // Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); // Rvs_SNP_Reads
        pChr = strtok(NULL,"\t"); // Total_SNP_Reads
        nSR = atoi(pChr); 

        bINS_Found = bDEL_Found = false;

        for (int i=nOffset-2; i < nOffset+3; i++)
        {
            pstKey->nOffset = i; pstKey->cAllele = acAllele[2];
            itr_INS = g_aINSList[nChromosome-1].find(pstKey);

            if (itr_INS!=g_aINSList[nChromosome-1].end())
            {                 
                if ((itr_INS->second == 1 && nSR > 11) || (nSR - itr_INS->second > 20)) 
                   continue;
                 
                nDiff_INS = abs(nOffset-i);nOffset_INS = pstKey->nOffset;                     
                bINS_Found=true; break;
            }
        }

        for (int i=nOffset-2; i < nOffset+3; i++)
        {
            pstKey->nOffset = i; pstKey->cAllele = acAllele[0];
            itr_DEL = g_aDELList[nChromosome-1].find(pstKey);

            if (itr_DEL!=g_aDELList[nChromosome-1].end())
            {   
                if ((itr_DEL->second == 1 && nSR > 11)||(nSR - itr_DEL->second > 20))
                     continue;

                nDiff_DEL = abs(nOffset-i); nOffset_DEL = pstKey->nOffset;
                bDEL_Found=true; break;
            }           
        }

        if (!bINS_Found && !bDEL_Found) continue;

        fprintf(g_pfOut,"%s\tS-%u:%s,%u\tI-",acChro,nOffset,acAllele,nSR);

        if (!bINS_Found) {fprintf(g_pfOut,"0\t");}            
        else{
            fprintf(g_pfOut,"%u:",nOffset_INS);
            fprintf(g_pfOut,"%c,%u",(*itr_INS->first).cAllele,itr_INS->second);
            fprintf(g_pfOut,"\t%u\t",nDiff_INS);
        }

        fprintf(g_pfOut,"D-");
        if (!bDEL_Found) {fprintf(g_pfOut,"0");}
        else{
            fprintf(g_pfOut,"%u:",nOffset_DEL);
            fprintf(g_pfOut,"%c,%u",(*itr_DEL->first).cAllele,itr_DEL->second);
            fprintf(g_pfOut,"\t%u\t",nDiff_DEL);
        }
        fprintf(g_pfOut,"\n"); 
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
    
    g_pfIn_SNP = fopen(argv[1],"r"); if (!g_pfIn_SNP) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto ExitRtn;}
    g_pfIn_INS = fopen(argv[2],"r"); if (!g_pfIn_INS) {fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto ExitRtn;}       
    g_pfIn_DEL = fopen(argv[3],"r"); if (!g_pfIn_DEL) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto ExitRtn;}
    g_pfOut = fopen(argv[4],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto ExitRtn;}

    GenerateList(g_pfIn_INS, g_aINSList);
    GenerateList(g_pfIn_DEL, g_aDELList); OutputData();

 ExitRtn:
    if (g_pfIn_SNP) fclose (g_pfIn_SNP); if (g_pfIn_INS) fclose (g_pfIn_INS); 
    if (g_pfIn_DEL) fclose (g_pfIn_DEL); if (g_pfOut) fclose (g_pfOut);

    return (EXIT_SUCCESS);

}
