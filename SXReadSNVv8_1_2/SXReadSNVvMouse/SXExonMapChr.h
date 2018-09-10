/* 
 * File:   SXExonMapChr.h
 * Author: Hwah Shih Yiew
 *
 * Created on November 26, 2009, 11:04 AM
 */

#ifndef _SXEXONMAPCHR_H
#define	_SXEXONMAPCHR_H

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include<vector>
#include <algorithm>
#include <string.h>
#define maxlen 1024

inline size_t Exon_rowParser(char * _source, const char * _delimiter, char ** _fields, size_t _max)
{
    size_t index = 0;
    char * pch = strtok(_source, _delimiter);
    while((pch != NULL) && (index < _max))
    {
        _fields[index++] = pch;
        pch = strtok(NULL, _delimiter);
    }
    return index;
}

struct Exon
{
    size_t start;
    size_t stop;
    char * exonId;

    inline Exon()
    {
        start = 0;
        stop = 0;
        exonId = NULL;
    }

    inline void insert(size_t _stop, char * _exonId, char * _ccds)
    {
        this->stop = _stop;

        /*
        unsigned sLen = strlen(_exonId) + strlen(_ccds);
        this->exonId = new char [sLen + 2];
        strcpy(this->exonId, _exonId);
        strcat(this->exonId, ",");
        strcat(this->exonId, _ccds);
        */

        unsigned sLen_exonId = strlen(_exonId);

        if (sLen_exonId > 1) 
        {
           unsigned sLen = sLen_exonId + strlen(_ccds);
           this->exonId = new char [sLen + 2];
           strcpy(this->exonId, _exonId);
           strcat(this->exonId, ",");
           strcat(this->exonId, _ccds);
        }
        else
        {          
           this->exonId = new char [strlen(_ccds)+ 1];
           strcpy(this->exonId, _ccds);
        }
    }

    ~Exon()
    {
        if (this->exonId != NULL)
            delete [] this->exonId;
    }
};


typedef struct _compareStartOffset{
    bool operator()(const Exon *ps1, const Exon *ps2)
    {
        //return ps1->start<ps2->start;

        if (ps1->start != ps2->start)
           return ps1->start < ps2->start;
        else
           return ps1->stop < ps2->stop;
    }
};

typedef std::vector<Exon*> EXON_VEC;
struct Exon_Chromosome
{
private:
    EXON_VEC exonMap;
    EXON_VEC::iterator exonIt;

public:
    inline void sortExon()
    {        
        sort(this->exonMap.begin(),this->exonMap.end(),_compareStartOffset());
    }

    inline void insertExon(size_t _start, size_t _stop, char * _exonId, char * _ccds)
    {
        Exon *pnewExon = new Exon();
        pnewExon->start = _start;
        pnewExon->insert(_stop, _exonId, _ccds);
        this->exonMap.push_back(pnewExon);        
    }

    inline size_t getSize()
    {
        return this->exonMap.size();
    }

    inline void printAll(const char* pcFile)
    {
        FILE* pf = fopen(pcFile,"w");
        for (this->exonIt = this->exonMap.begin(); this->exonIt != this->exonMap.end(); this->exonIt++)
            fprintf(pf, "%lu\t%lu\t%s\n", (*this->exonIt)->start, (*this->exonIt)->stop, (*this->exonIt)->exonId);
    }

    inline unsigned search(size_t _offset, char ** _exonId, unsigned _maxExon)
    {
        unsigned cnt = 0;
        
        EXON_VEC::iterator start_itr = this->exonMap.begin();

        if (start_itr == this->exonMap.end()) return cnt;
       
        unsigned int unTotalRecs = this->exonMap.size();
        unsigned int unPrevStartPos = 0;
        unsigned int unStartPos = unTotalRecs/2;
        unsigned int unLastPos = unTotalRecs;
        unsigned int unMaxPlus1 = unTotalRecs+1;
        int nLoop=0;

        do {
            start_itr = this->exonMap.begin()+unStartPos;

            if (_offset < (*start_itr)->start)
            {
                unLastPos = unStartPos; //-1;
                unStartPos = (unPrevStartPos+unLastPos)/2; //unStartPos /= 2;
            }
            else if (_offset > (*start_itr)->stop)
            {
                if (!(unLastPos >(unStartPos+1)))
                    break;

                unPrevStartPos = unStartPos;
                unStartPos=(unStartPos+unLastPos)/2;
            }
            else
            {
                break;
            }
            nLoop++;
        } while (unStartPos > 0 && unStartPos < unMaxPlus1 && nLoop < 100);
        

        if (!(nLoop < 100)){
            fprintf(stderr,"offset = %u\n", _offset);
            start_itr = this->exonMap.begin();
        }

        if (unStartPos > unTotalRecs) return cnt;
        if (unStartPos >(_maxExon-1)) start_itr -=_maxExon;//move _maxExon steps backwards to cater for offset repetition in exon lists
        else start_itr = this->exonMap.begin();
        
        for (this->exonIt=start_itr; this->exonIt != this->exonMap.end(); this->exonIt++)
        {
            if (_offset < (*this->exonIt)->start)
                break;

            if (!(_offset > (*this->exonIt)->stop))
            {
                _exonId[cnt++] = (*this->exonIt)->exonId;

                if (cnt == _maxExon)
                    break;
            }
        }
        return cnt;
    }
};


typedef struct _compareChrNum{
    bool operator()(unsigned s1, unsigned s2)
    {
        return s1<s2;
    }
};

typedef std::map<unsigned, Exon_Chromosome*, _compareChrNum> EXON_LIST;

bool GenerateExonMap(FILE * exonFile);
bool OutputExonID(unsigned int chrNum, unsigned int offset, FILE* pf);
void ClrExonMap();

#endif	/* _SXEXONMAPCHR_H */

