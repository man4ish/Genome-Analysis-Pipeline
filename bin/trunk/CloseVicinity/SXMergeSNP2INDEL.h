#ifndef _SXMERGESNP2INDEL_H
#define _SXMERGESNP2INDEL_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <map>

typedef struct _stKey{
    int nOffset;
}stKey;


typedef struct _stData
{
    int nFwdCnt;
    int nRvsCnt;
    int nReadCnt;
    unsigned int unReadDen;
    int nPEndCnt;
    float fAvgQryPos;
    float fAvgQScore;
    //char* pcdbSNP;

    _stData()
    {
        //pcdbSNP=NULL;
        nFwdCnt=nRvsCnt=nReadCnt=nPEndCnt=0;
        unReadDen=0;
        fAvgQryPos=fAvgQScore=0;
    }

    /*
    ~ _stData()
    {
        if (pcdbSNP) delete [] pcdbSNP;
    }
    */
}stData;


struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        return p1->nOffset < p2->nOffset;

    }
};


typedef std::map<stKey*,stData*,_compareKey> SNP_LIST;
typedef std::map<stKey*,int,_compareKey> GAP_LIST;


SNP_LIST g_aSNPList[24];
GAP_LIST g_aGapList[24];

#endif  /* _SXMERGESNP2INDEL_H */

