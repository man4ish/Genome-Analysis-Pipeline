/* 
 * File:   SSNpParser.h
 * Author: Hwah Shih Yiew
 *
 * Created on December 9, 2009, 5:35 PM
 */

#ifndef _SSNPPARSER_H
#define	_SSNPPARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <algorithm>
#include <map>

typedef struct _stKey{    
    unsigned int unGiNum;
    unsigned int unOffset;
    char cRefBase;
    char cVarBase;       
} stKey;

typedef struct _stSyno{
    char* pcGeneIDs;
    char acStatus[2];
    char acProtVar[2];

    _stSyno()
    {
        acProtVar[0]='\0';
    }
    
    ~_stSyno()
    {
        if (pcGeneIDs)
            delete[] pcGeneIDs;
    }

} stSyno;


typedef struct _compare{
    bool operator()(const stKey *p1, const stKey *p2)
    {
        if (p1->unGiNum != p2->unGiNum)
            return p1->unGiNum < p2->unGiNum;

        if (p1->unOffset != p2->unOffset)
            return p1->unOffset < p2->unOffset;

        if (p1->cRefBase != p2->cRefBase)
            return p1->cRefBase < p2->cRefBase;

        return p1->cVarBase < p2->cVarBase;
    }
};


typedef std::map<stKey*,stSyno*, _compare> SYNO_LIST;  //typedef std::vector<stSyno*> SYNO_LIST;


typedef struct _stStats
{
    unsigned int unTotalRec;
    unsigned int unS;
    unsigned int unNS;
    unsigned int unMM;
    unsigned int unNH;
    unsigned int unOR;

    _stStats()
    {
        unTotalRec=unS=unNS=unMM=unNH=unOR=0;
    }
} stStats;

inline char *strupr(char *s) { char *p = s; while (*p) { *p = toupper(*p); p++; } return s; }

inline size_t rowParser(char * _source, const char * _delimiter, char ** _fields, size_t _max)
{
    if (_max < 1)
        return 0;

    size_t index = 0;
    _fields[index++] = _source;

    while (*_source && index <= _max)
    {
        if (strchr(_delimiter, *_source))
        {
            // Found delimiter in current position
            *_source = '\0';
            if (strlen(_fields[index - 1]) == 0)
                index--;
            _fields[index++] = ++_source;
        }
        else
            _source++;
        if (*_source == '\n')
            *_source = '\0';
    }
    return index;
}


inline bool inRange(size_t _left, size_t _right, size_t _index)
{
    size_t left = _left;
    size_t right = _right;
    if (right < left)
    {
        size_t temp = left;
        left = right;
        right = temp;
    }
    if ((_index <= right) && (_index >= left))
        return true;
    return false;
}

inline double Percentage(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator*100;
}

void ProcessSSOutFile();
void PrintNewSNPList();
void PrintStatsRpt();
#endif	/* _SSNPPARSER_H */

