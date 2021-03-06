/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on June 1, 2010, 11:06 AM
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <map>

typedef struct _stKey
{
    int nStart;
    int nStop;
}stKey;

typedef struct _comparekey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
         if (p1->nStart != p2->nStart) return p1->nStart < p2->nStart;              
         return p1->nStop < p2->nStop;
    }
};


typedef struct _stData
{
   int nCnt;
}stData;


typedef std::map<stKey*, stData*, _comparekey> POLY_PDENS_LIST;

POLY_PDENS_LIST g_aPolyPDens[24];

FILE *g_pfIn=NULL, *g_pfOut=NULL, *g_pfPoly=NULL;// **g_apf=NULL; 
char g_cSample;

static void banner(char *argv[])
{
    printf("Synamatix SXAppendPDenMinus2Poly_SingleSampleMNV Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Poly Perfect Density filename > <Sample:- 1 or 2> <Output filename>\n\n");
}

void GenerateList()
{
    fprintf(stdout,"GenerateList\n");
    char acbuf[4096],acOut[4096]; char *pChr=NULL; stKey *pstKey=NULL; stData *pstData=NULL;
    int nChromosome; unsigned long ulTotalCnt;    
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfPoly))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateList-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfPoly)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;        

        strcpy(acOut,acbuf); pstKey = new stKey; pstData = new stData;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nStart = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstKey->nStop = atoi(pChr);
        pChr = strtok(NULL,"\t"); pstData->nCnt = atoi(pChr);

       if (pstData->nCnt < 1 ) { 
          fprintf(stdout,"%d\t%u\t%u\t%u\n",nChromosome+1,pstKey->nStart,pstKey->nStop,pstData->nCnt);
       }

        g_aPolyPDens[nChromosome-1][pstKey]=pstData;
    }
}


void AppendPerfectDen(char cChromosome, int nOffset, int nSNPCnt)
{
/*
    int nChroIdx=cChromosome-1; if (!g_apf[nChroIdx]) return;

    int nTotalPerfectRD=0; unsigned char acBuf[2];unsigned int unReads;
    
    nOffset--;
    fseek(g_apf[nChroIdx],nOffset*2,SEEK_SET);
    fread(acBuf,2,1,g_apf[nChroIdx]);

    unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];

    unReads = (unReads > 2)? unReads-2:0;  
    nTotalPerfectRD = unReads+nSNPCnt;

    fprintf(g_pfOut,"\t%d\n",nTotalPerfectRD);
*/
}


void AppendPerfectDenEx(char cChromosome, int nOffset,int nSNPCnt)
{ 
/* 
    int nChroIdx=cChromosome-1; if (!g_apf[nChroIdx]) return;

    int nTotalPerfectRD=0; unsigned char acBuf[2]; 
    unsigned int unRead=0, unTotalReads=0;

    nOffset--;

    //nOffset = nOffset-1;

    for (int i=0;i<2;i++)
    {
         fseek(g_apf[nChroIdx],(nOffset+i)*2,SEEK_SET);
         fread(acBuf,2,1,g_apf[nChroIdx]);

         unRead = acBuf[1];unRead<<=8;
         unRead|= acBuf[0];
         unTotalReads += unRead;
    }

    nTotalPerfectRD = (unsigned int)((double)((double)(unTotalReads)/2)+0.5);
    nTotalPerfectRD = (nTotalPerfectRD > 2)? nTotalPerfectRD-2:0;
    nTotalPerfectRD += nSNPCnt;

    fprintf(g_pfOut,"\t%d\n",nTotalPerfectRD);
*/
}


void AppendPolyPDens(char cChromosome, int nStart, int nStop)
{
    int nTotalPerfectRD = 0; int nChroIdx = cChromosome-1;
    stKey *pstKey = new stKey; pstKey->nStart = nStart; pstKey->nStop = nStop;

    POLY_PDENS_LIST::iterator itr=g_aPolyPDens[nChroIdx].find(pstKey);

    if (itr != g_aPolyPDens[nChroIdx].end()) nTotalPerfectRD=(*itr->second).nCnt;
        
    fprintf(g_pfOut,"\t%d\n",nTotalPerfectRD);     
}


void ProcessRecs()
{
    int nStart,nStop,nSNPCnt/*,nTotalReadDens*/,nNonSNPs, nReadStrength;
    
    char acbuf[1024],acOut[1024],*pChr=NULL, cChromosome;   
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {        
        if (++ulTest1 == 100000)
        {        
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }
 
        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;

        strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); nStart = atoi(pChr);         
        pChr = strtok(NULL,"\t"); nStop = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Nucleotide_Variant
        pChr = strtok(NULL,"\t");  //SNP_Density_1 
        pChr = strtok(NULL,"\t");  //SNP_Density_2
        pChr = strtok(NULL,"\t");  //Total_SNP_Density
        nSNPCnt = atoi(pChr);

/*
        pChr = strtok(NULL,"\t");  //Total_SNP_Reads
        nTotalSNPCnt = atoi(pChr);

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t");  //Total_Read_Density

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t");  //Total AvgQScore

        pChr = strtok(NULL,"\t");  //dbSnp
        pChr = strtok(NULL,"\t");  //ReadStrength   
        nReadStrength = atoi(pChr);

        bHomo = !(nReadStrength < g_fCutOff);

        if (!bHomo && (nTotalSNPCnt < g_nHetMinSR)) continue;
        if (bHomo && (nTotalSNPCnt < g_nHomMinSR)) continue;  
*/
        fprintf(g_pfOut,"%s",acOut);
     
        AppendPolyPDens(cChromosome,nStart,nStop);        
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
    if (g_pfPoly) fclose(g_pfPoly);  
    /* 
    iif (g_apf){
        for (int j=0;j<24; j++){
            if (g_apf[j]) fclose(g_apf[j]);
        }

        delete[] g_apf;
    }
    */
}


int main(int argc, char** argv) {

    if (argc < 5){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfPoly = fopen(argv[2],"r");
    if (!g_pfPoly){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_cSample = argv[3][0];
    
    g_pfOut = fopen(argv[4],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto Err_Exit;
    }
   
    GenerateList(); ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

