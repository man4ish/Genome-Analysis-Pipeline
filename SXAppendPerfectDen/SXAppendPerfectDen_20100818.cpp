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

FILE *g_pfIn=NULL, *g_pfOut=NULL, **g_apf=NULL; 


static void banner(char *argv[])
{
    printf("Synamatix SXAppendPerfectDen Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Perfect Density files path> <Output filename>\n\n");
}


void AppendPerfectDen(char cChromosome, int nOffset, int nSNPCnt)
{
    int nChroIdx=cChromosome-1; if (!g_apf[nChroIdx]) return;

    int nTotalPerfectRD=0; unsigned char acBuf[2];unsigned int unReads;
    
    nOffset--;
    fseek(g_apf[nChroIdx],nOffset*2,SEEK_SET);
    fread(acBuf,2,1,g_apf[nChroIdx]);

    unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];

    //unReads = (unReads > 2)? unReads-2:0;  
    nTotalPerfectRD = unReads+nSNPCnt;

    fprintf(g_pfOut,"\t%d\n",nTotalPerfectRD);
}


void ProcessRecs()
{
    int nOffset,nSNPCnt/*,nTotalReadDens*/,nNonSNPs, nReadStrength;
    
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

        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr);        
        pChr = strtok(NULL,"\t");  //Nucleotide_Variant
        pChr = strtok(NULL,"\t");  

        if (pChr[0] != '-')        //SNP_Density_1 
           nSNPCnt = atoi(pChr);
        else {
           pChr = strtok(NULL,"\t");  //SNP_Density_2
           nSNPCnt = atoi(pChr);
        }

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
     
        AppendPerfectDen(cChromosome,nOffset, nSNPCnt);
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

    if (argc < 3){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }
    
    g_apf = new FILE*[24];  char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_read_den_perfect",argv[2],j+1);
       g_apf[j] = fopen(acFile,"rb");

       if (!g_apf[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }
   
    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

