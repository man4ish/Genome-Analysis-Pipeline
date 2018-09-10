/*
 * File:   SXGeneMapChr.h
 * Author: jameswong
 * Modified/Enhanced by: syhwah
 * Created on October 21, 2009, 4:52 PM
 */

#ifndef _SXGENEMAPCHR_H
#define	_SXGENEMAPCHR_H

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <set>
#include <string.h>

inline size_t rowParser(char * _source, const char * _delimiter, char ** _fields, size_t _max)
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


struct ltstr
{
    inline bool operator()(const char * _lhs, const char * _rhs)
    {
        return (strcmp(_lhs, _rhs) < 0);
    }
};

struct Gene
{
    size_t start;
    size_t stop;
    char *geneId;

    inline Gene(size_t _start, size_t _stop)
    {
        start = _start;
        stop = _stop;
        geneId = NULL;        
    }

    inline void insert(char * _geneId)
    {
        this->geneId = _geneId;        
    }
};

struct GeneLtr
{
    inline bool operator()(const Gene * _lhs, const Gene * _rhs)
    {
        if (_lhs->start != _rhs->start)
            return (_lhs->start < _rhs->start);
        return (_lhs->stop < _rhs->stop);
    }
};

struct Chromosome
{
private:
    std::map<char*,char*, ltstr> geneIdSet;
    std::map<char*,char*, ltstr>::iterator geneIdIt;
    std::multiset<Gene*, GeneLtr> geneSet;
    std::multiset<Gene*, GeneLtr>::iterator geneIt;
    size_t maxGeneSize;

public:
    inline Chromosome()
    {
        maxGeneSize = 0;
    }

    inline void insertGene(size_t _start, size_t _stop, char *_geneId, char *_geneDesc)
    {
        this->geneIt = this->geneSet.insert(new Gene(_start, _stop));

        geneIdIt = geneIdSet.find(_geneId);
        if(geneIdIt == geneIdSet.end())
        {
            char *tmpGeneId = new char[strlen(_geneId)+1];
            strcpy(tmpGeneId, _geneId);
            char *tmpGeneDesc = new char[strlen(_geneDesc)+1];
            strcpy(tmpGeneDesc, _geneDesc);
            geneIdSet[tmpGeneId]=tmpGeneDesc;
            (*geneIt)->geneId = tmpGeneId;
        }
        else
            (*geneIt)->geneId = geneIdIt->first;

        size_t geneLen = _stop - _start + 1;
        if (geneLen > maxGeneSize)
            maxGeneSize = geneLen;
    }

    inline void printAll()
    {
        for (geneIt = geneSet.begin(); geneIt != geneSet.end(); geneIt++)
            fprintf(stderr, "%lu\t%lu\t%s\n", (*geneIt)->start, (*geneIt)->stop, (*geneIt)->geneId);
     
    }

    inline unsigned search(size_t _offset, char ** _geneId, char **_geneDesc, unsigned _maxGene)
    {
        std::set<char*, ltstr> fndGeneSet;
        unsigned cnt = 0;

        for (geneIt = geneSet.begin(); geneIt != geneSet.end(); geneIt++)
        {
            if (((*geneIt)->start <= _offset) && (*geneIt)->stop >= _offset)
                fndGeneSet.insert((*geneIt)->geneId);
            else if ((*geneIt)->start >= _offset)
                break;
        }

        std::set<char*, ltstr>::iterator fndGeneIt;
        for (fndGeneIt = fndGeneSet.begin(); fndGeneIt != fndGeneSet.end(); fndGeneIt++)
        {
            _geneId[cnt] = *fndGeneIt;

            this->geneIdIt = this->geneIdSet.find(*fndGeneIt);
            if (this->geneIdIt!= this->geneIdSet.end())
                _geneDesc[cnt]= this->geneIdIt->second;

            if (++cnt == _maxGene)
                break;
        }
        return cnt;
    }

    ~Chromosome()
    {
        for (geneIt = geneSet.begin(); geneIt != geneSet.end(); geneIt++)
            delete [] *(this->geneIt);

        for (this->geneIdIt = this->geneIdSet.begin(); this->geneIdIt != this->geneIdSet.end(); this->geneIdIt++)
        {         
            delete[] this->geneIdIt->first; delete[] this->geneIdIt->second;
        }        
    }
};



bool GenerateGeneMap(FILE *geneFile);
bool OutputGeneIDnDesc(unsigned int chrNum, unsigned int offset, FILE* pf);
void ClrGeneMap();
#endif	/* _SXGENEMAPCHR_H */

