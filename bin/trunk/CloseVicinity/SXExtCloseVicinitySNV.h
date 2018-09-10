#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string.h>
#include <ctype.h>

typedef struct _stKey
{
    int nOffset;
    char cAllele;
}stKey;

/*
typedef struct _stData
{
    int nCnt;
    char acAllele[4];
    int anSReads[4];

    _stData()
    {
        nCnt =0;
        memset(acAllele,0,4);
        for (int i=0; i<4;i++){ anSReads[i]=0; }   
    }
}stData;
*/


struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        if (p1->nOffset != p2->nOffset) 
    	   return p1->nOffset < p2->nOffset;
        
        return p1->cAllele < p2->cAllele;
    }
};


//typedef std::map<stKey*, stData*, _compareKey> SNV_LIST;

typedef std::map<stKey*,int, _compareKey> SNV_LIST;

SNV_LIST g_aINSList[24], g_aDELList[24];
FILE *g_pfIn_SNP, *g_pfIn_INS, *g_pfIn_DEL, *g_pfOut; 


