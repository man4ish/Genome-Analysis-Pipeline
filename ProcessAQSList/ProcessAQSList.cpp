
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include<map>

typedef struct _stKey{
    char cChromosome;
    int nOffset;
    char acAllele[4];
} stKey;


typedef struct _stData{

}stData;

typedef struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        if (p1->cChromosome != p2->cChromosome)
            return p1->cChromosome < p2->cChromosome;

        if (p1->nOffset != p2->nOffset)
            return p1->nOffset < p2->nOffset;

        if (p1->acAllele[0] != p2->acAllele[0])
            return p1->acAllele[0] < p2->acAllele[0];

        return p1->acAllele[2] < p2->acAllele[2]; //strcmp(p1->acAllele,p2->acAllele);
    }
};


typedef std::map<stKey*,char,_compareKey> SXSAMPLE_LIST;

SXSAMPLE_LIST g_SXSampleList; 
FILE *g_pfInput1=NULL, *g_pfInput2=NULL, *g_pfOutput=NULL;


void GenerateKeyList()
{
    char acbuf[10240]; char *pChr=NULL;

    stKey *pstKey=NULL;  SXSAMPLE_LIST::iterator itr=NULL;

    while(!feof(g_pfInput1))
    {
        if (!fgets(acbuf,sizeof(acbuf),g_pfInput1)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C')
           continue;

        pstKey = new stKey;

        pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0]))
            pstKey->cChromosome = atoi(pChr);
        else if (pChr[0]=='X')
            pstKey->cChromosome = 23;
        else if(pChr[0]=='Y')
            pstKey->cChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);

           g_SXSampleList[pstKey]=1;
    }

}

void OutputRecList()
{
    char acbuf[10240],acOut[10240]; char *pChr=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
    stKey *pstKey=new stKey;  SXSAMPLE_LIST::iterator itr=NULL;

    while(!feof(g_pfInput2))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfInput2)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;
        if (acbuf[0] == 'C') fprintf(g_pfOutput,"%s\n",acbuf);

        strcpy(acOut,acbuf); pChr = strtok(acbuf,"\t");

        if (!isalpha(pChr[0])) pstKey->cChromosome = atoi(pChr);
        else if (pChr[0]=='X') pstKey->cChromosome = 23;
        else if(pChr[0]=='Y') pstKey->cChromosome = 24;

        pChr = strtok(NULL,"\t"); pstKey->nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); strcpy(pstKey->acAllele,pChr);

        itr = g_SXSampleList.find(pstKey);

        if (itr != g_SXSampleList.end()) fprintf(g_pfOutput,"%s\n",acOut);          
        
    }

}


int main(int argc, char *argv[] )
{
     g_pfInput1 = fopen(argv[1],"r"); if (!g_pfInput1){printf("Failed to open %s ...\n",argv[1]); exit(9);}
     g_pfInput2 = fopen(argv[2],"r"); if (!g_pfInput2){printf("Failed to open %s ...\n",argv[2]); exit(9);}
     g_pfOutput = fopen(argv[3],"w"); if (!g_pfOutput){printf("Failed to open %s ...\n",argv[3]); exit(9);}

     GenerateKeyList(); OutputRecList(); 
}
