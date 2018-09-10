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

FILE *g_pfIn=NULL, *g_pfOut=NULL,**g_apfMaskedDen=NULL,**g_apfReadDen=NULL,**g_apfReptDen=NULL; 
char g_cSNVType; float g_fCutOff,g_fMeans;
int g_nHetMinSR, g_nHomMinSR;

static void banner(char *argv[])
{
    printf("Synamatix SXPrintZygosityCNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Means> <Read Density files path> <Repeat Density files path> "
           "<Masked Density files path> <Output filename>\n\n");
}


void OutputLocalCNV(char cChromosome, int nOffset)
{
    int nChroIdx=cChromosome-1; if (!g_apfMaskedDen[nChroIdx]) return;

    unsigned int unReads, unTotalRD=0, unTotalCnt=0, unReadDenVals=0, unReptDenVals=0; 
    unsigned char acBuf[2];

    //char acMask[5000]={'\0'}, acRead[5000]={'\0'}, acRept[5000]={'\0'}; 

    int nFranking = 100;

    bool bPrint = false; if (nOffset == 122872) bPrint = true;

    nOffset--;
    nOffset -= nFranking;

    for (int nIdx=0; nIdx< ((nFranking*2)+1); nIdx++)
    {        
        fseek(g_apfMaskedDen[nChroIdx],nOffset*2,SEEK_SET);
        fread(acBuf,2,1,g_apfMaskedDen[nChroIdx]);

        unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
        //sprintf(acMask+strlen(acMask),"%u,",unReads); 
        
        if (unReads == 0){
            unTotalCnt++;

            fseek(g_apfReadDen[nChroIdx],nOffset*2,SEEK_SET);
            fread(acBuf,2,1,g_apfReadDen[nChroIdx]);
            unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
            unTotalRD += unReads; unReadDenVals += unReads;          
            
            //if (bPrint) fprintf(stdout,"%d\t%u\n",nOffset*2,unReads);
 
            //sprintf(acRead+strlen(acRead),"%u;",unReads);  

            fseek(g_apfReptDen[nChroIdx],nOffset*2,SEEK_SET);
            fread(acBuf,2,1,g_apfReptDen[nChroIdx]);
            unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
            unTotalRD += unReads; unReptDenVals += unReads;

            if (bPrint) fprintf(stdout,"%d\t%u\n",nOffset*2,unReads); 
            //sprintf(acRept+strlen(acRept),"%u,",unReads);
        }
    //    else
    //    {
    //        sprintf(acRead+strlen(acRead),"0,"); sprintf(acRept+strlen(acRept),"0,");
    //    }
 
        nOffset++; 
    }
      
    fprintf(g_pfOut,"\t%.2f\t%u\t%u\t%u\n",((double)unTotalRD/(double)unTotalCnt)/(g_fMeans/2),unTotalCnt,unReadDenVals,unReptDenVals);   
    //fprintf(g_pfOut,"\t%s\t%s\t%s\t%u\t%u\t%.2f\t%u\t%u\n",
    //                acMask, acRead, acRept, unTotalRD, unTotalCnt,
    //                ((double)unTotalRD/(double)unTotalCnt)/(g_fMeans/2),unReadDenVals,unReptDenVals);
}


void ProcessRecs()
{
    int nOffset,nTotalSNPCnt/*,nTotalReadDens*/,nNonSNPs, nReadStrength;
    
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

        strcpy(acOut,acbuf);

        if (acbuf[0] == 'C'){
           fprintf(g_pfOut,"%s\tMappable_Bases\tMappable_Read_Densities\tMappable_Repeat_Densities\n",acOut);
           continue;                 
        }

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        nOffset = atoi(pChr);
/*
        pChr = strtok(NULL,"\t");  //Nucleotide_Variant
        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
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
     
        OutputLocalCNV(cChromosome,nOffset);
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

    if (g_apfMaskedDen){
        for (int j=0;j<24; j++){
            if (g_apfMaskedDen[j]) fclose(g_apfMaskedDen[j]);
        }

        delete[] g_apfMaskedDen;
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
    /* 
    g_cSNVType = toupper(argv[4][0]);
    if (g_cSNVType != 'S' && g_cSNVType != 'I' && g_cSNVType != 'D'){
        fprintf(stdout, "Invalid SNVType Specified %c ...\n",g_cSNVType); goto Err_Exit;
    }
    */
    g_apfReadDen = new FILE*[24]; char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_read_den",argv[3],j+1);
       g_apfReadDen[j] = fopen(acFile,"rb");

       if (!g_apfReadDen[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }

    g_apfReptDen = new FILE*[24]; 

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_rept_den",argv[4],j+1);
       g_apfReptDen[j] = fopen(acFile,"rb");

       if (!g_apfReptDen[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }

    g_apfMaskedDen = new FILE*[24]; 

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_mask_den",argv[5],j+1);
       g_apfMaskedDen[j] = fopen(acFile,"rb");

       if (!g_apfMaskedDen[j]){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }
    }


    g_pfOut = fopen(argv[6],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[6]); goto Err_Exit;
    }

    //g_fCutOff = atof(argv[2]), g_nHetMinSR = atoi(argv[3]), g_nHomMinSR = atoi(argv[4]), 
    g_fMeans = atof(argv[2]);
    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

