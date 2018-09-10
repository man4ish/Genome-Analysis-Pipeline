#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL,**g_apfReadDen=NULL,**g_apfReptDen=NULL,*g_pfOut=NULL; 
float g_fMean=0.00,g_fXCoverage=0.00,g_fYCoverage=0.00;
float fPEScore,fSNPReadScore,fStrandScore,fAvgQS,fRS,fReadDenScore,fSNPScore;


static void banner(char *argv[])
{
    printf("Synamatix SXAppendScoreSNP_SingleSample Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Read Density files path> "
           "<Repeat Density files path> <Output filename>\n\n");
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


void OutputConfidentScore(char cChromosome, int nOffset, int nFwdSNPCnt, 
                          int nRvsSNPCnt, int nSNPCnt, int nTotalRD, int nPECnt)
{
    int nChroIdx=cChromosome-1;

    unsigned int unReads, unTotalRD=0, unTotalCnt=0, unReadDenVals=0, unReptDenVals=0;
    unsigned int unMBPReads=0, unReptReads=0, unMBPCnt=0, unMBPVals=0;
    unsigned char acBuf[2];

    int nFranking = 100;

    nOffset--;
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
    
    fprintf(g_pfOut,"-\t-\t-\t-\t%u\t%u\t%u\t%u\t%u\t",unTotalCnt,unReadDenVals,unReptDenVals,unMBPCnt,unMBPVals);

    int nMappableBases=unTotalCnt,nMappableReptDens=unReptDenVals,nReptDenScore,nSNPReadScore,nStrandScore;

    fPEScore = (nPECnt > 2)?100:(float)nPECnt*(100.00/3.00);
    nSNPReadScore = (nSNPCnt > 4)?100:nSNPCnt*20;
    nStrandScore = ((nFwdSNPCnt > 1) && (nRvsSNPCnt > 1))?100:0;
    fRS = (float)nSNPCnt/(float)nTotalRD*100.00;
    fReadDenScore = (fRS > 10.00)?100.00:fRS*10.00;

    if (nMappableReptDens==0)
        nReptDenScore =20;
    else
        nReptDenScore = ((float)nMappableReptDens/(float)nMappableBases < 66.00)?20:0; //16+(0.5*100) => 66.00

    fSNPScore = (fPEScore+(float)nSNPReadScore+(float)nStrandScore+fReadDenScore+(float)nReptDenScore)/420.00;

    fprintf(g_pfOut,"%.2f\t%d\t%d\t%.2f\t%d\t%.4f\n",
                     fPEScore,nSNPReadScore,nStrandScore,fReadDenScore,nReptDenScore,fSNPScore);

 
}


void ProcessRecs()
{
    int nOffset,nTotalSNPCnt,nTotalRD,nNonSNPs, nReadStrength;
    int nFwdSNPCnt,nRvsSNPCnt,nSNPCnt,nPECnt;

    char acbuf[40960],acOut[40960],*pChr=NULL, cChromosome;
    bool bHomo;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;
        if (acbuf[0] == 'C') {
           fprintf(g_pfOut,"%s\tMappable_Bases\tMappable_Read_Densities\tMappable_Repeat_Densities"
                           "\tMappable_Bases_Repeat\tRepeat_Densities_MBP\tPEScore\tSNPReadScore"
                           "\tStrandScore\tReadDenScore\tReptDenScore\tConfidentScore\n",acbuf); 
           continue;
        }  

        strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t"); //Chromosome

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); //GiNumber
        pChr = strtok(NULL,"\t"); //Offset
        nOffset = atoi(pChr);

        pChr = strtok(NULL,"\t"); //Nucleotide_Variant
        pChr = strtok(NULL,"\t");  nFwdSNPCnt = atoi(pChr); //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  nRvsSNPCnt = atoi(pChr); //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  nSNPCnt = atoi(pChr);    //SNP_Reads                 
        pChr = strtok(NULL,"\t");  nTotalRD = atoi(pChr);   //Total_Read_Density
        pChr = strtok(NULL,"\t");  nPECnt = atoi(pChr);     //PEnd_Count 

        fprintf(g_pfOut,"%s\t",acOut);
        OutputConfidentScore(cChromosome,nOffset,nFwdSNPCnt,nRvsSNPCnt,nSNPCnt,nTotalRD,nPECnt);
   }
}


int main(int argc, char** argv) { 

    if (argc < 5){
        fprintf(stdout,"Parameters not sufficient...\n\n"); goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    //g_fMean = atof(argv[2]); g_fXCoverage=g_fMean/2.00; g_fYCoverage=g_fMean;

    g_apfReadDen = new FILE*[24];  char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_read_den",argv[2],j+1);
       g_apfReadDen[j] = fopen(acFile,"rb");

       if (!g_apfReadDen[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }


    g_apfReptDen = new FILE*[24];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_rept_den",argv[3],j+1);
       g_apfReptDen[j] = fopen(acFile,"rb");

       if (!g_apfReptDen[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    } 

    g_pfOut = fopen(argv[4],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

