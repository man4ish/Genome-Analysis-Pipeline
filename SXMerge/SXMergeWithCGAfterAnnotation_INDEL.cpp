/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on May 25, 2010, 10:30 AM
 */

#include "SXMergeWithCGAfterAnnotation_INDEL.h"

/*
 * 
 */

SNP_LIST g_SNPList;
FILE *g_pfFile1, *g_pfFile2, *g_pfCommon, *g_pfUnique1, *g_pfUnique2;


void GenerateSNPList()
{
    fprintf(stdout,"GenerateSNPList\n");
    char acbuf[40960],acOut[40960]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    SNP_LIST::iterator itr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfFile1))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateSNPList-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFile1)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
           continue;

        strcpy(acOut,acbuf); 

        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) pstKey->cChromosome = atoi(pChr);
        else if (pChr[0]=='X') pstKey->cChromosome = 23;
        else if(pChr[0]=='Y') pstKey->cChromosome = 24;

        pChr = strtok(NULL,"\t"); 
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);
        
        pstData = new stData;
        pstData->pcData = new char[strlen(acOut)+1];
        strcpy(pstData->pcData,acOut);

        g_SNPList[pstKey]=pstData;   
    }
}


void AnalyseSNPData()
{
    fprintf(stdout,"AnalyseSNPData\n");   
    char acbuf[40960], acOut[40960]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    SNP_LIST::iterator itr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;  

    while(!feof(g_pfFile2))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"AnalyseSNPData-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFile2)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0] == '\t'){
           continue;
        }

        strcpy(acOut,acbuf); 

        pstKey = new stKey;

        pChr = strtok(acbuf,"\t"); //locus
        pChr = strtok(NULL,"\t");  //ploidy
        pChr = strtok(NULL,"\t");  //haplotype
        pChr = strtok(NULL,"\t");  //chromosome

        if (!isalpha(pChr[0])) pstKey->cChromosome = atoi(pChr);
        else if (pChr[0]=='X') pstKey->cChromosome = 23;
        else if(pChr[0]=='Y') pstKey->cChromosome = 24;
        /* 
        pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);
        */

        pChr = strtok(NULL,"\t");   //begin

        //stKey.uOffset = atoi(pChr);

        pChr = strtok(NULL,"\t");  //end
        pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); //varType
        //pChr = strtok(NULL,"\t"); //reference
        //pstKey->acAllele[0] = pChr[0];
        pChr = strtok(NULL,"\t"); //alleleSeq
        pstKey->acAllele[0] = pChr[0];         

        itr = g_SNPList.find(pstKey);

        if (itr==g_SNPList.end()){
            fprintf(g_pfUnique2,"%s\n",acOut);
        }
        else{
            (*itr->second).bUnique = false;
            (*itr->second).pcCGData = new char[strlen(acOut)+1];
            strcpy((*itr->second).pcCGData,acOut);
        }        

        delete pstKey;
    }
}


void OutputSNPData()
{
    AnalyseSNPData();

    fprintf(stdout,"OutputSNPData\n"); 
    int nTotalSNPCnt, nTotalReadDensity;      
    FILE *pf=NULL; SNP_LIST::iterator itr = g_SNPList.begin();
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;  

    while(itr!=g_SNPList.end())
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        pf = ((*itr->second).bUnique)? g_pfUnique1 : g_pfCommon;

        fprintf(pf,"%s",(*itr->second).pcData);
        
        if ((*itr->second).pcCGData) fprintf(pf,"\t%s\n",(*itr->second).pcCGData);
        else fprintf(pf,"\n"); 

        itr++; 
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

    GenerateSNPList(); OutputSNPData();

 ExitRtn:
    if (g_pfFile1) fclose (g_pfFile1); 
    if (g_pfFile2) fclose (g_pfFile2);
    if (g_pfCommon) fclose (g_pfCommon);
    if (g_pfUnique1) fclose (g_pfUnique1); 
    if (g_pfUnique2) fclose (g_pfUnique2);

    g_SNPList.clear();

    return (EXIT_SUCCESS);
}

