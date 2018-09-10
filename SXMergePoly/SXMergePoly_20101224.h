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
    char cChromosome;
    int nOffset;
    char acAllele[4];
}stKey;


typedef struct _stData
{
    int nSNPCnt;
    int nReadDensity;
    float fAvgQScore;
    char *pcdbSNP;
    bool bUnique;
    int nSNPCnt2;
    int nReadDensity2;
    float fAvgQScore2;

    _stData()
    {
        pcdbSNP=NULL;
        bUnique=true;
        nSNPCnt2=nReadDensity2=0;
        fAvgQScore2=0;
    }

    ~ _stData()
    {
        if (pcdbSNP) delete [] pcdbSNP;
    }
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

        return p1->acAllele[2] < p2->acAllele[2];
    }
};


typedef std::map<stKey*,stData*,_compareKey> SNP_LIST;

#endif	/* _SXMERGE_H */

