/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on June 1, 2010, 11:06 AM
 */

#include "SXPrintCNV.h"

/*
 * 
 */


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
    int nChroIdx=cChromosome-1; 

    unsigned int unReads, unTotalRD=0, unTotalCnt=0, unReadDenVals=0, unReptDenVals=0; 
    unsigned int unMBPReads=0, unReptReads=0, unMBPCnt=0, unMBPVals=0;
    unsigned char acBuf[2];

    int nFranking = 100;

    nOffset--;
    nOffset -= nFranking;

    for (int nIdx=0; nIdx< ((nFranking*2)+1); nIdx++)
    {        
        memcpy(acBuf,&(g_stMaskedDen[nChroIdx].pcData[nOffset*2]),2);
        unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
        
        if (unReads == 0){
            unTotalCnt++;

            memcpy(acBuf,&(g_stReadDen[nChroIdx].pcData[nOffset*2]),2);
            unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
            unTotalRD += unReads; unReadDenVals += unReads; unMBPReads=unReads;         

            memcpy(acBuf,&(g_stReptDen[nChroIdx].pcData[nOffset*2]),2);
            unReads = acBuf[1];unReads<<=8; unReads|= acBuf[0];
            unTotalRD += unReads; unReptDenVals += unReads; unReptReads=unReads;

            if (unMBPReads ==0 && unReptReads!=0) {
               unMBPCnt++; unMBPVals +=unReads;
            }
        }
 
        nOffset++; 
    }
      
    fprintf(g_pfOut,"\t%.2f\t%u\t%u\t%u\t%u\t%u\n",
            ((double)unTotalRD/(double)unTotalCnt)/(g_fMeans/(double)2),
            unTotalCnt,unReadDenVals,unReptDenVals,unMBPCnt,unMBPVals);   
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
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u M\n", ulTest2);
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

    //int rc;
    for (int j=0;j<24; j++){
        if (g_stReadDen[j].fd > 0){
            munmap(g_stReadDen[j].pcData,g_stReadDen[j].FileSize);
            close(g_stReadDen[j].fd);  
        }
    }

    for (int j=0;j<24; j++){
        if (g_stReptDen[j].fd > 0){
            munmap(g_stReptDen[j].pcData,g_stReptDen[j].FileSize);
            close(g_stReptDen[j].fd);
        }
    }

    for (int j=0;j<24; j++){
        if (g_stMaskedDen[j].fd > 0){
            munmap(g_stMaskedDen[j].pcData,g_stMaskedDen[j].FileSize);
            close(g_stMaskedDen[j].fd);
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

    struct stat st; char acFile[1024]; 

    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_read_den",argv[3],j+1);
        
       g_stReadDen[j].fd=open(acFile, O_RDONLY, S_IRUSR);

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
       sprintf(acFile,"%s/Chr%02d_rept_den",argv[4],j+1);

       g_stReptDen[j].fd=open(acFile, O_RDONLY, S_IRUSR); 

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


    for (int j=0;j<24; j++)
    {
       sprintf(acFile,"%s/Chr%02d_mask_den",argv[5],j+1);

       g_stMaskedDen[j].fd=open(acFile, O_RDONLY, S_IRUSR); 

       if (g_stMaskedDen[j].fd == -1){ 
           fprintf(stdout,"Failed to open %s ...",acFile); goto Err_Exit;
       }

       fstat(g_stMaskedDen[j].fd,&st); g_stMaskedDen[j].pcData = new char[st.st_size];
       g_stMaskedDen[j].pcData=(char*)mmap(NULL,st.st_size,PROT_READ,MAP_PRIVATE,g_stMaskedDen[j].fd,0);

       if (g_stMaskedDen[j].pcData==(void*)-1){
           fprintf(stderr, "Error: %s mmap failed: %s\n",acFile, strerror(errno));goto Err_Exit;
       }
       g_stMaskedDen[j].FileSize=st.st_size;
    }


    g_pfOut = fopen(argv[6],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[6]); goto Err_Exit;
    }

    
    g_fMeans = atof(argv[2]);
    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

