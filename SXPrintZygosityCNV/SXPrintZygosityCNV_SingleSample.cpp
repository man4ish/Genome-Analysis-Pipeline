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

static void banner(char *argv[])
{
    printf("Synamatix SXPrintZygosityCNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <CutOff> <Output filename>\n\n");

    //printf("\t<Input filename> <CutOff> <Type:- S,I,D> <Output filename>\n\n");
}


void OutputLocalCNV(char cChromosome, int nOffset)
{
    int nChroIdx=cChromosome-1; if (!g_apf[nChroIdx]) return;

    int nTotalRD=0; unsigned char acBuf[2];unsigned int unReads;

    int nFranking = 100;    
    nOffset--;
    nOffset -= nFranking; 

    for (int nIdx=0; nIdx< ((nFranking*2)+1); nIdx++)
    {
        fseek(g_apf[nChroIdx],nOffset*2,SEEK_SET);
        fread(acBuf,2,1,g_apf[nChroIdx]);

        unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
        nTotalRD += unReads;
        nOffset++;
    }

    fprintf(g_pfOut,"\t%.2f\n",(((double)(nTotalRD))/201)/g_fMeans);
}


void ProcessRecs()
{
    int nOffset,nTotalSNPCnt,nTotalReadDens, nNonSNPs;
    char acbuf[40960],acOut[40960],*pChr=NULL, cChromosome;
    bool bHomo;

    while(!feof(g_pfIn))
    {
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
        pChr = strtok(NULL,"\t");  //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  //Rvs_SNP_Reads

        pChr = strtok(NULL,"\t");  //Total_SNP_Reads
        nTotalSNPCnt = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Total_Read_Density
        nTotalReadDens = atoi(pChr);

        /*  
        if (g_cSNVType == 'I'){
            nNonSNPs = nTotalReadDens - nTotalSNPCnt;
            if (nNonSNPs > nTotalSNPCnt) nTotalSNPCnt = nNonSNPs;
        }
        */

        bHomo=((int)((double)(nTotalSNPCnt*100)/(double)(nTotalReadDens)))>g_fCutOff;

        fprintf(g_pfOut,"%s\t%s\n",acOut,bHomo?"hom":"het");
        //OutputLocalCNV(cChromosome,nOffset);
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

    if (argc < 3){ //4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    /*  
    g_cSNVType = toupper(argv[3][0]);
    if (g_cSNVType != 'S' && g_cSNVType != 'I' && g_cSNVType != 'D'){
        fprintf(stdout, "Invalid SNVType Specified %c ...\n",g_cSNVType); goto Err_Exit;
    }
    */ 

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    g_fCutOff = atof(argv[2]), //g_fMeans = atof(argv[3]);
    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

