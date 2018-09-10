#include <stdio.h>
#include <stdlib.h>
#include <map>

typedef struct _stKey
{
    int nOffset;
}stKey;


typedef struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
    	return p1->nOffset < p2->nOffset;
    }
};


typedef std::map<stKey*, int, _compareKey> SNV_LIST;

SNV_LIST g_aINSList[24], g_aDELList[24];
FILE *g_pfIn_SNP, *g_pfIn_INS, *g_pfIn_DEL, *g_pfOut; 


