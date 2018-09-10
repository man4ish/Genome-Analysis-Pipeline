#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <map>

typedef struct _stKey
{
     char cChromo;
     int nOffset;
} stKey;


typedef struct _compare{
    bool operator()(const stKey *pKey1, const stKey *pKey2)
    {
        if (pKey1->cChromo != pKey2->cChromo)
           return pKey1->cChromo < pKey2->cChromo;    

        return pKey1->nOffset < pKey2->nOffset;  
    }
};


typedef std::map<stKey*, short, _compare> SNV_List;


FILE *g_pfIn=NULL,*g_pfFilterIn=NULL,*g_pfOut=NULL,*g_pfFilterOut=NULL; 

SNV_List g_FilterList; SNV_List::iterator g_Itr;

static void banner(char *argv[])
{
    printf("Synamatix SXFilterSNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Input Filtered filename> <Output filename> < Output Filtered filename>\n\n");
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfFilterIn) fclose(g_pfFilterIn);
    if (g_pfOut) fclose(g_pfOut); if (g_pfFilterOut) fclose(g_pfFilterOut);
}


void GenerateFilteredList()
{
    char acbuf[1024],*pChr=NULL; stKey *pstKey=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfFilterIn))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"GenerateFilteredList - Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFilterIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pstKey = new stKey; 
        
        pChr = strtok(acbuf,"\t"); pstKey->cChromo = atoi(pChr);
        pChr = strtok(NULL,"\t");  pstKey->nOffset = atoi(pChr);        
        
        g_Itr = g_FilterList.find(pstKey);

        if (g_Itr == g_FilterList.end()) g_FilterList[pstKey]=1;        
        else delete[] pstKey;                 
   } 
}


void ProcessRecs()
{
    char acbuf[1024],acOut[1024],*pChr=NULL; stKey *pstKey=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"ProcessRecs - Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pstKey = new stKey; strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t"); pstKey->cChromo = atoi(pChr);
        pChr = strtok(NULL,"\t");  pstKey->nOffset = atoi(pChr);

        g_Itr = g_FilterList.find(pstKey);

        if (g_Itr == g_FilterList.end()) fprintf(g_pfOut,"%s\n",acOut);        
        else fprintf(g_pfFilterOut,"%s\n",acOut);

        delete[] pstKey;        
   }    
}


int main(int argc, char** argv) {

    if (argc < 5){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfFilterIn = fopen(argv[2],"r");
    if (!g_pfFilterIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    g_pfFilterOut = fopen(argv[4],"w");
    if (!g_pfFilterOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto Err_Exit;
    }

    GenerateFilteredList(); ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

