#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL, *g_pfOut=NULL;

typedef struct _stRecIn{
    char cChromo;
    int nOffset;
    char acData[1024];
    bool bToBeExtracted;

    _stRecIn()
    {
       cChromo=0; nOffset=0; bToBeExtracted = false;  
    }  

}stRecIn;


static void banner(char *argv[])
{
    printf("Synamatix SXExtractConseq Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Output filename>\n\n");
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


void ProcessRecs()
{
    char acbuf[1024],acOut[1024],*pChr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; 

    stRecIn stPrevRec, stCurrRec;

    if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){return;} 

    strcpy(acOut,acbuf);
    pChr = strtok(acbuf,"\t"); stPrevRec.cChromo = pChr[0];
      
    pChr = strtok(NULL,"\t");  stPrevRec.nOffset = atoi(pChr);
    strcpy(stPrevRec.acData,acOut);      

    while(!feof(g_pfIn))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %lu x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;
        strcpy(stCurrRec.acData,acbuf);

        pChr = strtok(acbuf,"\t"); stCurrRec.cChromo = pChr[0];
        pChr = strtok(NULL,"\t");  stCurrRec.nOffset = atoi(pChr);
        
        if (stCurrRec.cChromo==stPrevRec.cChromo && stCurrRec.nOffset==stPrevRec.nOffset+1)
        {
           fprintf(g_pfOut,"%s\n",stPrevRec.acData); 
           stPrevRec.bToBeExtracted = true;
        }     
        else 
        {
           if (stPrevRec.bToBeExtracted) fprintf(g_pfOut,"%s\n",stPrevRec.acData);   
           stPrevRec.bToBeExtracted = false;            
        }   
        stPrevRec.nOffset = stCurrRec.nOffset; strcpy(stPrevRec.acData, stCurrRec.acData);          
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

    g_pfOut = fopen(argv[2],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

