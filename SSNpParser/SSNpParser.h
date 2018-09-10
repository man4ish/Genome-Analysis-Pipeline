/* 
 * File:   SSNpParser.h
 * Author: Manish
 *
 * Created on December 9, 2009, 5:35 PM
 */

#ifndef _SSNPPARSER_H
#define	_SSNPPARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <algorithm>
#include <map>

typedef struct _stKey{    
    unsigned int unGiNum;
    unsigned int unOffset;
    char cRefBase;
    char cVarBase;      
   
   _stKey()
   {
   	unGiNum=unOffset=0;
        cRefBase=cVarBase=0;     
   } 
} stKey;

typedef struct _stSyno{
    char* pcGeneIDs;
    char acStatus[2];
    char acProtVar[2];
    //bool bNonStop;

    _stSyno()
    {
        acProtVar[0]='\0';
        // = true;
    }
    
    ~_stSyno()
    {
        if (pcGeneIDs)
            delete[] pcGeneIDs;
    }

} stSyno;


struct _compare{
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


inline double Division(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator;
}

char g_acMonth[12][4]= {"Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec"};

inline void PrintRptDateTime(FILE *pf)
{
    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    fprintf(pf,"%02d-%s-20%02d; %02d:%02d:%02d\n",
               ptm->tm_mday,g_acMonth[ptm->tm_mon],ptm->tm_year-100,
               ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
}


inline void PrintRemarksParam(FILE *pfParam)
{
    FILE *pf = fopen("Param.ini","r");
    if (!pf) return;

    char acbuf[1024]; size_t unSize=sizeof(acbuf);
    char *pChr=NULL;

    while(!feof(pf))
    {
        if (!fgets(acbuf,unSize,pf)){break;}
        if (acbuf[0]=='#') continue;

        pChr = strtok(acbuf,"=");

        if (strcmp(pChr,"FRemarks")==0)
        {
            pChr = strtok(NULL,"\n");

            if (pChr)
                fprintf(pfParam,"%s",pChr);
        }
    }

    fprintf(pfParam,"\n");
    if (pf) {fclose(pf);}
}

void ProcessSSOutFile();
void PrintNewSNPList();
void PrintStatsRpt();
//void OutputType(FILE *pf, char cRefBase, char cVarBase);
//bool IsNonStop(char *pcBases, char* pcSbjSeq,int nStart);
#endif	/* _SSNPPARSER_H */

