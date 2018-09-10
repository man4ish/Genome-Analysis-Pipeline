/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on June 1, 2010, 11:06 AM
 */

#include "SXPrintCNVNS.h"


/*
 * 
 */


static void banner(char *argv[])
{
    printf("Synamatix SXPrintCNVMNV Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <SampleID> <Means> <Read Density files path> <Repeat Density files path> "
           "<Output filename>\n\n");
}


void OutputLocalCNV(char cChromosome, int nStart, int nStop)
{
    int nChroIdx=cChromosome-1;

    unsigned int unReads, unTotalRD=0, unTotalCnt=0, unReadDenVals=0, unReptDenVals=0; 
    unsigned int unMBPReads=0, unReptReads=0, unMBPCnt=0, unMBPVals=0;
    unsigned char acBuf[2];

    int nOffset,nFranking = 100;

    for (int i=nStart; i < (nStop+1); i++) 
    {
         nOffset = i-1;
         nOffset -= nFranking;

         for (int nIdx=0; nIdx< ((nFranking*2)+1); nIdx++)
         {         
              memcpy(acBuf,&(g_stReadDen[nChroIdx].pcData[nOffset*2]),2);            
              unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
              unTotalRD += unReads; unReadDenVals += unReads; unMBPReads=unReads;         
            
              memcpy(acBuf,&(g_stReptDen[nChroIdx].pcData[nOffset*2]),2);
              unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
              unTotalRD += unReads; unReptDenVals += unReads; unReptReads=unReads;

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
    int nStart, nStop;
    
    char acbuf[1024],acOut[1024],*pChr=NULL, cChromosome;
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

        strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
        else if (pChr[0]=='X') cChromosome = 23;
        else if(pChr[0]=='Y') cChromosome = 24;
        else cChromosome = 0;

        pChr = strtok(NULL,"\t"); nStart = atoi(pChr);
        pChr = strtok(NULL,"\t"); nStop = atoi(pChr);
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

        OutputLocalCNV(cChromosome,nStart,nStop);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);

    for (int j=0;j<24; j++){
        if (g_stReadDen[j].fd > 0) {
            munmap(g_stReadDen[j].pcData,g_stReadDen[j].FileSize);
            close(g_stReadDen[j].fd);
        }
    }
    
    
    for (int j=0;j<24; j++){
        if (g_stReptDen[j].fd > 0) {
           munmap(g_stReptDen[j].pcData,g_stReptDen[j].FileSize);
           close(g_stReptDen[j].fd);
        }
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

    struct stat st; char acFile[1024];

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/%s.seq%d.read_den_all",argv[4],g_pcSampleID,j+1);
       g_stReadDen[j].fd = open(acFile,O_RDONLY,S_IRUSR);

       if (g_stReadDen[j].fd == -1){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }

       fstat(g_stReadDen[j].fd,&st); g_stReadDen[j].pcData = new char[st.st_size];
       g_stReadDen[j].pcData=(char*)mmap(NULL,st.st_size,PROT_READ,MAP_PRIVATE,g_stReadDen[j].fd,0);

       if (g_stReadDen[j].pcData==(void*)-1){
          fprintf(stderr, "Error: %s mmap failed: %s\n",acFile, strerror(errno));goto Err_Exit;
       }
       g_stReadDen[j].FileSize=st.st_size;
    }


    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/%s.seq%d.rept_den",argv[5],g_pcSampleID,j+1);
       g_stReptDen[j].fd = open(acFile,O_RDONLY,S_IRUSR);

       if (g_stReptDen[j].fd == -1){
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }

       fstat(g_stReptDen[j].fd,&st); g_stReptDen[j].pcData = new char[st.st_size];
       g_stReptDen[j].pcData=(char*)mmap(NULL,st.st_size,PROT_READ,MAP_PRIVATE,g_stReptDen[j].fd,0);

      if (g_stReptDen[j].pcData==(void*)-1){
         fprintf(stderr, "Error: %s mmap failed: %s\n",acFile, strerror(errno));goto Err_Exit;
      } 
      g_stReptDen[j].FileSize=st.st_size;
    }

    g_pfOut = fopen(argv[6],"w");

    if (!g_pfOut){ fprintf(stdout,"Failed to open %s ...\n",argv[5]); goto Err_Exit; }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

