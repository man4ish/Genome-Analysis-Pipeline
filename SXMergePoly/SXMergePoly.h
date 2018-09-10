/* 
 * File:   SXMerge.h
 * Author: Hwah Shih Yiew
 *
 * Created on May 25, 2010, 10:31 AM
 */

#ifndef _SXMERGE_H
#define	_SXMERGE_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <map>

typedef struct _stKey{
    int nOffset;
    char acAllele[4];

   _stKey()
   {
      memset(acAllele,0,4);
   } 
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
    char* pcdbSNP;   

    _stData()
    {
        pcdbSNP=NULL;
        nFwdCnt=nRvsCnt=nReadCnt=nPEndCnt=0;
        unReadDen=0;
        fAvgQryPos=fAvgQScore=0; 
    }

    ~ _stData()
    {
        if (pcdbSNP) delete [] pcdbSNP;
    }
}stData;


struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        if (p1->nOffset != p2->nOffset)
            return p1->nOffset < p2->nOffset;

        if (p1->acAllele[0] != p2->acAllele[0])
            return p1->acAllele[0] < p2->acAllele[0];

        return p1->acAllele[2] < p2->acAllele[2];
    }
};


typedef std::map<stKey*,stData*,_compareKey> SNP_LIST;

SNP_LIST g_aSNPList[24];

#endif	/* _SXMERGE_H */

