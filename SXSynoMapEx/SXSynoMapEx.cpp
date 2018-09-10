/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on March 02, 2011, 4:01 AM
 */

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
/*
 * 
 */


FILE **g_apfHsRef;
char **g_apcHsRef;
int g_anHsRefSize[24];
unsigned int *g_punGINums=NULL;


unsigned int g_aunGINums_G36[24] = {89161185,89161199,89161205,89161207,51511721,89161210,
                                   89161213,51511724,89161216,89161187,51511727,89161190,
                                   51511729,51511730,51511731,51511732,51511734,51511735,
                                   42406306,51511747,51511750,89161203,89161218,89161220};


unsigned int g_aunGINums_G37[24] = {224589800,224589811,224589815,224589816,224589817,224589818,
                                   224589819,224589820,224589821,224589801,224589802,224589803,
                                   224589804,224589805,224589806,224589807,224589808,224589809,
                                   224589810,224589812,224589813,224589814,224589822,224589823};


static void banner(char *argv[])
{
    printf("Synamatix SXSynoMapEx Copyright 2009 Synamatix Sdn Bhd\n");
    printf("Build %s %s\n",__DATE__,__TIME__);
    printf("Usage: %s <Human Ref. files path> <Genome Ver.> ><input file name> <output file name>\n\n", argv[0]);
}


void CloseRefFiles()
{
    for (short i=0;i<24;i++) {
        if (g_apfHsRef[i]){ fclose(g_apfHsRef[i]); g_apfHsRef[i]=NULL;}
    }
}


void LoadHRefSeq()
{
    struct stat st; size_t file_size;
    char acbuf[130]; char *pChr=NULL;
    int nLen,nTotalRecs; unsigned int unIdx;

    for (short i=0;i<24;i++)
    {
        fstat(fileno(g_apfHsRef[i]), &st);
        file_size = (size_t)st.st_size;
        fgets(acbuf, sizeof(acbuf),g_apfHsRef[i]);

        nLen = strlen(acbuf); nTotalRecs = file_size - nLen;

        g_apcHsRef[i] = new char[nTotalRecs]; unIdx=0;

        while(!feof(g_apfHsRef[i]))
        {
            if (!fgets(acbuf, sizeof(acbuf),g_apfHsRef[i])) break;

            pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
            pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

            nLen = strlen(acbuf); memcpy(g_apcHsRef[i]+unIdx,acbuf,nLen); unIdx+=nLen;
        }
        g_apcHsRef[i][unIdx]=0; g_anHsRefSize[i]=unIdx;

       fprintf(stdout,"Human Ref.File %d => ",i+1);
       for (int j=0; j < 70; j++){ fprintf(stdout,"%c",g_apcHsRef[i][j]); }
       fprintf(stdout,"...\n");
   }

   CloseRefFiles();
}


bool OpenHRefFiles(const char *pcRefFilePath)
{
    char acFileName[1024];

    for (short i=0;i<24;i++)
    {
        if (i < 22){
            sprintf(acFileName,"%s/hs_ref_GRCh37_chr%d.fa",pcRefFilePath,i+1);
         }
        else if(i == 22)
            sprintf(acFileName,"%s/hs_ref_GRCh37_chrX.fa",pcRefFilePath);
        else
            sprintf(acFileName,"%s/hs_ref_GRCh37_chrY.fa",pcRefFilePath);

        g_apfHsRef[i] = fopen(acFileName,"r");

        if (g_apfHsRef[i]){
           fprintf(stdout,"Touch %s ...\n",acFileName);
        }
        else{
            fprintf(stdout,"Failed to open %s ....\n",acFileName);
            CloseRefFiles(); return false;
        }
    }
    return true;
}


void Initialization()
{
    g_apfHsRef = new FILE*[24]; g_apcHsRef = new char*[24];

    for (short i=0;i<24;i++){
        g_apfHsRef[i]=NULL; g_apcHsRef[i]=NULL;
    }
}


int main(int argc, char** argv) {

    Initialization();
      
    if (argc <5) {printf("Parameters required...\n");banner(argv); exit(9);}     

    if (!OpenHRefFiles(argv[1])) exit(9);

    g_punGINums = (strcmp(argv[2],"37.1")==0)?g_aunGINums_G37:g_aunGINums_G36;
  
    FILE *pfSyn = fopen(argv[3],"r"); 
    if (!pfSyn){fprintf(stderr,"Failed to open %s ...\n",argv[3]); exit(9);}

    FILE *pf=fopen(argv[4],"w");

    if (!pf){fprintf(stderr,"Failed to open %s ...\n",argv[4]); exit(9);}

    fprintf(stdout,"Load Human Reference Sequence ...\n");

    LoadHRefSeq();  
 
    char acbuf[256]; unsigned int nBufSize = sizeof(acbuf); 
    char *pChr; char *lpc[5]; int nIdx;  char acSeq[90];
    int nOffset;

    while(!feof(pfSyn))
    {
        if (!fgets(acbuf,nBufSize,pfSyn)) break;
        
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        nIdx=0; pChr = strtok(acbuf,"|");

        while (pChr)
        {
             lpc[nIdx++]= pChr;
             if (nIdx > 4) break;  pChr = strtok(NULL,"|");
        }
             
        for (int i=0; i <24; i++){              
             if ((unsigned int)atoi(lpc[1]) == g_punGINums[i]){
                nOffset = atoi(lpc[2])-45;                                
                memcpy(&acSeq[0],&g_apcHsRef[i][nOffset],89); 
                acSeq[44]=lpc[3][1]; acSeq[89]=0;                 
                fprintf(pf,">%s|%s|%s|%s|%s\n%s\n",lpc[0],lpc[1],lpc[2],lpc[3],lpc[4],acSeq);
                break;
             } 
        }

        memset(lpc,0,sizeof(lpc));        
    }

    if (pfSyn) fclose(pfSyn); if (pf) fclose(pf);

    return (EXIT_SUCCESS);
}

