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
/*
 * 
 */

FILE *g_pfIn=NULL, *g_pfOut=NULL, **g_apfReadDen=NULL,**g_apfReptDen=NULL;
float g_fMeans; const char *g_pcSampleID;

static void banner(char *argv[])
{
    printf("Synamatix SXPrintZygosityCNV_SingleSampleMNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <SampleID> <Means> <Read Density files path> "
           "<Repeat Density files path> <Output filename>\n\n");
}


void OutputLocalCNV(char cChromosome, int nStart, int nStop)
{
    int nChroIdx=cChromosome-1; 

    unsigned int unReads, unTotalRD=0, unTotalCnt=0, unReadDenVals=0, unReptDenVals=0;
    unsigned int unMBPReads=0, unReptReads=0, unMBPCnt=0, unMBPVals=0;
    unsigned char acBuf[2];

    int nOffset,nFranking = 30;

    for (int i=nStart; i < (nStop+1); i++)
    {
         nOffset=i-1;
         nOffset -= nFranking;

         for (int nIdx=0; nIdx< ((nFranking*2)+1); nIdx++)
         {
              fseek(g_apfReadDen[nChroIdx],nOffset*2,SEEK_SET);
              fread(acBuf,2,1,g_apfReadDen[nChroIdx]);
              unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
              unTotalRD += unReads; unReadDenVals += unReads; unMBPReads=unReads;

              fseek(g_apfReptDen[nChroIdx],nOffset*2,SEEK_SET);
              fread(acBuf,2,1,g_apfReptDen[nChroIdx]);
              unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0]; unReptReads=unReads;
              unTotalRD += unReads; unReptDenVals += unReads;
         
              if (unMBPReads ==0 && unReptReads!=0) {
                  unMBPCnt++; unMBPVals +=unReads;
              }

              nOffset++;
         }
    }

    int nBps = (nStop - nStart)+1;     
    unReadDenVals = (unsigned int)(((double)unReadDenVals/(double)nBps)+0.5);
    unReptDenVals = (unsigned int)(((double)unReptDenVals/(double)nBps)+0.5);
    unMBPCnt = (unsigned int)(((double)unMBPCnt/(double)nBps)+0.5);
    unMBPVals = (unsigned int)(((double)unMBPVals/(double)nBps)+0.5); 

    fprintf(g_pfOut,"\t-\t%u\t%u\t%u\t%u\t%u\n",unTotalCnt,unReadDenVals,unReptDenVals,unMBPCnt,unMBPVals);
}


void ProcessRecs()
{
    int nStart,nStop,nTotalSNPCnt,nTotalReadDens;
    char acbuf[1024],acOut[1024],*pChr=NULL, cChromosome;
    float fRS;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %lu x 0.1M\n", ulTest2);
            ulTest1=0;
        }
 
        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;

        strcpy(acOut,acbuf); fprintf(g_pfOut,"%s",acOut);
        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); nStart = atoi(pChr);
        pChr = strtok(NULL,"\t"); nStop = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Nucleotide_Variant

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t"); //Fwd_SNP_Reads,Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  //SNP_Reads
        nTotalSNPCnt = atoi(pChr);

        pChr = strtok(NULL,"\t");  //_Read_Density
        nTotalReadDens = atoi(pChr);

        fRS=((float)nTotalSNPCnt)/((float)nTotalReadDens);
          
        if (fRS > 0.30 && fRS < 0.76) fprintf(g_pfOut,"\thet");
        else if (!(fRS < 0.76)) fprintf(g_pfOut,"\thom");
        else fprintf(g_pfOut,"\t-");
 
        OutputLocalCNV(cChromosome,nStart,nStop);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);

    if (g_apfReadDen){
        for (int j=0;j<24; j++){
            if (g_apfReadDen[j]) fclose(g_apfReadDen[j]);
        }

        delete[] g_apfReadDen;
    }

    if (g_apfReptDen){
        for (int j=0;j<24; j++){
            if (g_apfReptDen[j]) fclose(g_apfReptDen[j]);
        }

        delete[] g_apfReptDen;
    }
}


int main(int argc, char** argv) {

    if (argc < 7){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }
       
    g_pcSampleID = argv[2]; g_fMeans = atof(argv[3]);

    g_apfReadDen = new FILE*[24];  char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/%s.seq%d.read_den_all",argv[4],g_pcSampleID,j+1);
       g_apfReadDen[j] = fopen(acFile,"rb");

       if (!g_apfReadDen[j]){
           fprintf(stdout,"Failed to open %s ...\n",acFile); goto Err_Exit;
       }
    }

    g_apfReptDen = new FILE*[24];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/%s.seq%d.rept_den",argv[5],g_pcSampleID,j+1);
       g_apfReptDen[j] = fopen(acFile,"rb");

       if (!g_apfReptDen[j]){
           fprintf(stdout,"Failed to open %s ...\n",acFile); goto Err_Exit;
       }
    }

    g_pfOut = fopen(argv[6],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[6]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

