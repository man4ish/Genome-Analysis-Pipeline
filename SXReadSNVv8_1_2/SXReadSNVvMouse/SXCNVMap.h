/* 
 * File:   SXCNVMap.h
 * Author: Hwah Shih Yiew
 *
 * Created on February 2, 2010, 5:13 PM
 */

#ifndef _SXCNVMAP_H
#define	_SXCNVMAP_H

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <set>
#include <string.h>

inline size_t CNV_rowParser(char * _source, const char * _delimiter, char ** _fields, size_t _max)
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

typedef struct _stCNV
{
    size_t start;
    size_t stop;

    inline _stCNV(size_t _start, size_t _stop)
    {
        start = _start;
        stop = _stop;
    }
}stCNV;


struct CNVLtr
{
    inline bool operator()(const stCNV * _lhs, const stCNV * _rhs)
    {
        if (_lhs->start != _rhs->start)
            return (_lhs->start < _rhs->start);
        return (_lhs->stop < _rhs->stop);
    }
};

struct CNV_Chromosome
{
private:
    //std::map<char*,char*, ltstr> CNVIdSet;
    //std::map<char*,char*, ltstr>::iterator CNVIdIt;
    std::multiset<stCNV*, CNVLtr> CNVSet;
    std::multiset<stCNV*, CNVLtr>::iterator CNVIt;

public:

    inline void insertCNV(size_t _start, size_t _stop)
    {
        /*this->CNVIt = */this->CNVSet.insert(new stCNV(_start, _stop));
    }

    inline void printAll(const char* pcFile)
    {
        FILE* pf = fopen(pcFile,"w");

        //int m=CNVSet.size();
        for (CNVIt = CNVSet.begin(); CNVIt != CNVSet.end(); CNVIt++)
            fprintf(pf, "%lu\t%lu\n", (*this->CNVIt)->start, (*this->CNVIt)->stop);

        fclose(pf);
    }

    inline bool search(size_t _offset)
    {
        for (CNVIt = CNVSet.begin(); CNVIt != CNVSet.end(); CNVIt++)
        {
            if (((*CNVIt)->start <= _offset) && (*CNVIt)->stop >= _offset)
                return true;
            else if ((*CNVIt)->start > _offset)
                return false;
        }
        return false;
    }

    ~CNV_Chromosome()
    {
        for (CNVIt = CNVSet.begin(); CNVIt != CNVSet.end(); CNVIt++)
            delete [] *(this->CNVIt);
    }
};

bool GenerateCNVMap(FILE *geneFile);
void OutputCNVSearch(unsigned int chrNum, unsigned int offset, FILE* pf);
void ClrCNVMap();

#endif	/* _SXCNVMAP_H */

