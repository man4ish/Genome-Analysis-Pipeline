#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string.h>
#include <ctype.h>

typedef struct _stKey
{
    int nOffset;
}stKey;


struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        return p1->nOffset < p2->nOffset;
    }
};


typedef std::map<stKey*, int, _compareKey> SNV_LIST;

SNV_LIST g_aSNVList[24];
FILE *g_pfInPoly, *g_pfIn, *g_pfOut, *g_pfOutFiltered;


void GeneratePolybaseList()
{
    fprintf(stdout,"GeneratePolybaseList\n");
    char acbuf[4096],acInfo[1024],*pChr=NULL;
    stKey *pstKey=NULL;
    int nChromosome = 0;
    SNV_LIST::iterator itr;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfInPoly))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"GenerateList-Total Recs Processed = %lu M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInPoly)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
        {
           continue;
        }

        //strcpy(acOut,acbuf); 
        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        pChr = strtok(NULL,"\t"); 
        strcpy(acInfo,pChr); pChr = strtok(acInfo,":");

        pstKey->nOffset = atoi(pChr+2);       

        itr = g_aSNVList[nChromosome-1].find(pstKey);

        if (itr == g_aSNVList[nChromosome-1].end())
           g_aSNVList[nChromosome-1][pstKey]=1;
    }
}


void FilterData()
{
    char acbuf[4096],acOut[4096], *pChr=NULL;
    int nChromosome = 0;
    stKey *pstKey=new stKey; SNV_LIST::iterator itr;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"FilterData-Total Recs Processed = %lu M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0]==0||acbuf[0] == '#' || acbuf[0] == 'C') continue;

        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome = 23;
        else if(pChr[0]=='Y') nChromosome = 24;

        //pChr = strtok(NULL,"\t");// giNumber
        pChr = strtok(NULL,"\t");
        pstKey->nOffset = atoi(pChr);

        itr = g_aSNVList[nChromosome-1].find(pstKey);

        if (itr == g_aSNVList[nChromosome-1].end())
           fprintf(g_pfOut,"%s\n",acOut);
        else   
           fprintf(g_pfOutFiltered,"%s\n",acOut);
    }
}


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input Polybase filename> <Input filename-SNP List> <Output filename-new SNP List> <Output filename-SNP data filtered>\n\n");
       exit(9);
    }

    g_pfInPoly = fopen(argv[1],"r"); if (!g_pfInPoly) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto ExitRtn;}
    g_pfIn = fopen(argv[2],"r"); if (!g_pfIn) {fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto ExitRtn;}
    g_pfOut = fopen(argv[3],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto ExitRtn;}
    g_pfOutFiltered = fopen(argv[4],"w"); if (!g_pfOutFiltered) {fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto ExitRtn;}

    GeneratePolybaseList(); FilterData();

ExitRtn:

    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut); if (g_pfOutFiltered) fclose(g_pfOutFiltered);
}
