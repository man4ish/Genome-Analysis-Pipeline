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

FILE *g_pfIn=NULL, *g_pfOut=NULL;
int g_nXCoverage;

static void banner(char *argv[])
{
    printf("Synamatix SXPZygoINDEL Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Coverage> <Output filename>\n\n");
}


void ProcessRecs()
{
    int nOffset,nTotalSNPCnt,nTotalPerfectReadDens;
    float fRS;
    char acbuf[1024],*pChr=NULL, *pcLast=NULL, cChromosome;    

    while(!feof(g_pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;

        //strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t"); fprintf(g_pfOut,"%s",pChr); 

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); fprintf(g_pfOut,"\t%s",pChr);
        nOffset = atoi(pChr);

        pChr = strtok(NULL,"\t"); fprintf(g_pfOut,"\t%s",pChr); //Nucleotide_Variant
        pChr = strtok(NULL,"\t"); fprintf(g_pfOut,"\t%s",pChr); //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); fprintf(g_pfOut,"\t%s",pChr); //Rvs_SNP_Reads 
        pChr = strtok(NULL,"\t"); fprintf(g_pfOut,"\t%s",pChr); nTotalSNPCnt = atoi(pChr); 
        pChr = strtok(NULL,"\t"); pcLast=pChr;
    
        while (pChr)
        {
            pcLast=pChr;
            pChr = strtok(NULL,"\t");
            if (pChr) fprintf(g_pfOut,"\t%s",pcLast);
        }  
          
        nTotalPerfectReadDens = atoi(pcLast)-g_nXCoverage;

        if (nTotalSNPCnt > nTotalPerfectReadDens) nTotalPerfectReadDens = nTotalSNPCnt;

        fRS=((float)nTotalSNPCnt)/((float)nTotalPerfectReadDens);

        if (fRS > 0.30 && fRS < 0.76) fprintf(g_pfOut,"\thet\n");
        else if (!(fRS < 0.76)) fprintf(g_pfOut,"\thom\n");
        else fprintf(g_pfOut,"\t-\n"); 

    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv) {

    if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    if (atoi(argv[2]) > 15) g_nXCoverage = (int)((atof(argv[2])/10.00)+0.9);
    else g_nXCoverage = 0;     

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

