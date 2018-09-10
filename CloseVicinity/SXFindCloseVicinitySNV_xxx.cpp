#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "SXFindCloseVicinitySNV.h"


typedef struct _stKey
{
    int nOffset;
}stKey;


typedef struct _stData
{
    char *pcOut;

    ~_stData()
    {
         if (pcOut) delete[] pcOut; 
    }     
}stData;


typedef struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        return p1->nOffset < p2->nOffset;
    }
};


typedef std::map<stKey*,int,_compareKey> SNVTAG_List;

SNVTAG_LIST g_aTagList[24]; SNVTAG_LIST::iterator g_itr;

FILE *g_pfIn, *g_pfOut;


void GenerateList()
{
    char acbuf[1024],acOut[1024],*pChr=NULL; 
    stKey *pstKey; stData *pstData;
    int nChromosome;

    while (!feof(pfIn))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr='\0';
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr='\0';

        if (acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0]==0) continue;

        pstData = new stData;
        pstData->pcOut = new char[strlen(acData)+1];
        strcpy(pstData->pcOut,acData);

        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);  

        g_itr = g_aTagList[nChromosome-1].find(pstKey);

        if (g_itr == g_aTagList[nChromosome-1].end()){
            g_aTagList[nChromosome-1][pstKey] = pcData;
        }  
        else{
            delete [] pstData->pcOut;
        }
    }
}


void OutputData()
{     
    for (int nChroIdx=0; nChroIdx<24; nChroIdx++ )
    {
        for (g_itr = g_aTagList[nChroIdx].begin(); g_itr != g_aTagList[nChroIdx].end(); g_itr++)
        {
            if (srstr((*g_itr->second).pcOut),"SNPS")
            {
                for (int i=(*g_itr->first).nOffset-3; i < (*g_itr->first).nOffset+4; i++) 
                {
                }
            }
        }
    }
}


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        fprintf(stdout,"Parameters not sufficient...\n");

        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename> <Output filename>\n");

        exit(9);
    }

    unsigned long ulTest1=0, ulTest2=0;
    FILE *pfIn, *pfOut; 

    pfIn = fopen(argv[1],"r"); if (!pfIn) { fprintf(stdout,"Failed to opend %s ...\n", argv[1]); goto ExitRtn;}
    pfOut = fopen(argv[2],"w"); if (!pfOut) {fprintf(stdout,"Failed to opend %s ...\n", argv[3]); goto ExitRtn;}

    GenerateList(); OutputData();

ExitRtn:
    if (pfIn) fclose(pfIn); if (pfOut) fclose(pfOut);
}
