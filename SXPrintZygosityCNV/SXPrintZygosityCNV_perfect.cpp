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

FILE *g_pfIn=NULL, *g_pfOut=NULL, **g_apf=NULL; char g_cSNVType; float g_fCutOff,g_fMeans;
int g_nHetMinSR, g_nHomMinSR;

static void banner(char *argv[])
{
    printf("Synamatix SXPrintZygosityCNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <CutOff> <Min. Supporting Reads for het.> <Min. Supporting Reads for hom.> "
           "<Means> <Read Density files path> <Output filename>\n\n");
}


void OutputLocalCNV(char cChromosome, int nOffset)
{
    int nChroIdx=cChromosome-1; if (!g_apf[nChroIdx]) return;

    int nTotalRD=0; unsigned char acBuf[2];unsigned int unReads;

    nOffset -= 100;

    for (int nIdx=0; nIdx<201; nIdx++)
    {
        nOffset--;
        fseek(g_apf[nChroIdx],nOffset*2,SEEK_SET);
        fread(acBuf,2,1,g_apf[nChroIdx]);

        unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
        nTotalRD += unReads;
    }

    fprintf(g_pfOut,"\t%.2f\n",(((double)(nTotalRD))/201)/g_fMeans);
}


void ProcessRecs()
{
    int nOffset,nTotalSNPCnt/*,nTotalReadDens*/,nNonSNPs, nReadStrength;
    char acbuf[1024],acOut[1024],*pChr=NULL, cChromosome;
    bool bHomo;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u x 1.0M\n", ulTest2);
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

        pChr = strtok(NULL,"\t");
        nOffset = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Nucleotide_Variant
        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t");  //Total_SNP_Reads
        nTotalSNPCnt = atoi(pChr);

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t");  //Total_Read_Density
        //nTotalReadDens = atoi(pChr);

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t");  //Total AvgQScore

        pChr = strtok(NULL,"\t");  //dbSnp 
        pChr = strtok(NULL,"\t");  //ReadStrength   
        nReadStrength = atoi(pChr);

        //bHomo=((int)((double)(nTotalSNPCnt*100)/(double)(nTotalReadDens)))>g_fCutOff;

        bHomo = !(nReadStrength < g_fCutOff);

        if (!bHomo && (nTotalSNPCnt < g_nHetMinSR)) continue;
        if (bHomo && (nTotalSNPCnt < g_nHomMinSR)) continue;  

        fprintf(g_pfOut,"%s\t%s",acOut,bHomo?"hom":"het");
        OutputLocalCNV(cChromosome,nOffset);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);

    if (g_apf){
        for (int j=0;j<24; j++){
            if (g_apf[j]) fclose(g_apf[j]);
        }

        delete[] g_apf;
    }
}


int main(int argc, char** argv) {

    if (argc < 8){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }
    /* 
    g_cSNVType = toupper(argv[4][0]);
    if (g_cSNVType != 'S' && g_cSNVType != 'I' && g_cSNVType != 'D'){
        fprintf(stdout, "Invalid SNVType Specified %c ...\n",g_cSNVType); goto Err_Exit;
    }
    */
    g_apf = new FILE*[24];  char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_read_den_perfect",argv[6],j+1);
       g_apf[j] = fopen(acFile,"rb");

       if (!g_apf[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }

    g_pfOut = fopen(argv[7],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[7]); goto Err_Exit;
    }

    g_fCutOff = atof(argv[2]), g_nHetMinSR = atoi(argv[3]), g_nHomMinSR = atoi(argv[4]), g_fMeans = atof(argv[5]);
    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

