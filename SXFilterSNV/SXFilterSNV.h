#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <map>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>


typedef struct _stKey
{
     char nChromo;
     int nOffset;
} stKey;


struct _compare{
    bool operator()(const stKey *pKey1, const stKey *pKey2)
    {
        if (pKey1->nChromo != pKey2->nChromo)
           return pKey1->nChromo < pKey2->nChromo;

        return pKey1->nOffset < pKey2->nOffset;
    }
};


typedef std::map<stKey*, short, _compare> SNV_List;


#pragma pack(1)
typedef struct _stSNV
{
    int nChromosome;
    unsigned char acOffset[4];
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand; //0-rvs, 1-fwd
} stSNV;

#pragma pack()


int g_fd;

FILE *g_pfIn=NULL,*g_pfFilterIn=NULL,*g_pfOut=NULL,*g_pfFilterOut=NULL;

SNV_List g_FilterList; SNV_List::iterator g_Itr;

