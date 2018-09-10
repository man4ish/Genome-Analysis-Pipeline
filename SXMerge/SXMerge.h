#ifndef _SXMERGEV2_H
#define _SXMERGEV2_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define LIST_SIZE 1000000

typedef struct _stData{
    int nOffset;
    char acAllele[4];
    int nSNPCnt;
    int nReadDensity;
    float fAvgQScore;
    char *pcdbSNP;
    int nFwdSR;
    int nRvsSR;
    int nPEndCnt;
    float fAvgQryPos;
    bool bUnique;
    int nSNPCnt2;
    int nReadDensity2;
    float fAvgQScore2;
    int nFwdSR2;
    int nRvsSR2;
    int nPEndCnt2;
    float fAvgQryPos2;  

    _stData()
    {
        pcdbSNP=NULL;
        bUnique=true;
        nSNPCnt2=nReadDensity2=nFwdSR2=nRvsSR2=nPEndCnt2=0;
        fAvgQScore2=fAvgQryPos2=0;
    }

    ~ _stData()
    {
        if (pcdbSNP) delete [] pcdbSNP;
    }

}stData;


typedef struct _stRecList{
    stData *pData;
    int nUsed;
    int nTotal;    

    _stRecList()
    {
        pData = (stData*)malloc(sizeof(stData)*LIST_SIZE);    
        nUsed = 0;
        nTotal = LIST_SIZE;
    }

    ~_stRecList()
    {
       if (pData) free(pData);    
    }  
}stRecList;


#endif
