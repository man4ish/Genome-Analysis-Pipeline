/* 
 * File:   SXAnnotate.h
 * Author: Hwah Shih Yiew
 *
 * Created on June 2, 2010, 2:36 PM
 */

#ifndef _SXANNOTATE_H
#define	_SXANNOTATE_H

#include <string.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/mman.h>
#include <map>

using namespace std;

#define REC_SIZE 4096


typedef struct _stBinnedData
{
    int nOffset;
    int nIdx;
} stBinnedData;

typedef struct _stBinnedIdx
{
    stBinnedData *pBinnedData;
    int nUsedSubscript;
    int nTotalSubscript;

    _stBinnedIdx()
    {
        pBinnedData = (stBinnedData*)malloc(sizeof(stBinnedData)*REC_SIZE);
        nUsedSubscript=0;
        nTotalSubscript = REC_SIZE;
    }

    ~ _stBinnedIdx()
    {
        free(pBinnedData);
    }
} stBinnedIdx;


typedef struct _stExceptionData
{
    unsigned int nStart;
    unsigned int nStop;
    char *pcData;

    _stExceptionData()
    {
        nStart=nStop=0;
        pcData=NULL;
    }

    ~_stExceptionData()
    {
        if (pcData) delete [] pcData;
    }
}stExceptionData;


typedef struct _stExceptionIdx
{
    stExceptionData *pExceptionData;
    int nUsedSubscript;
    int nTotalSubscript;

    _stExceptionIdx()
    {
        pExceptionData = (stExceptionData*)malloc(sizeof(stExceptionData)*REC_SIZE);
        nUsedSubscript=0;
        nTotalSubscript = REC_SIZE;
    }

    ~ _stExceptionIdx()
    {
        free(pExceptionData);
    }
} stExceptionIdx;


typedef struct _stAnnotateRecs
{
    char *pcData;
    int nTotalSize;
    int nCurrSize;    

    _stAnnotateRecs()
    {
        pcData = (char*)malloc(REC_SIZE);
        nCurrSize=0;
        nTotalSize = REC_SIZE;
    }

    ~_stAnnotateRecs()
    {
        if (pcData) free(pcData);
    }
} stAnnotateRecs;


class SXAnnotate
{
public:
    SXAnnotate();
    ~SXAnnotate();

    char* GetRecord(char cChar, size_t qstart, size_t qstop);
    bool LoadDataSource(const char* pcbinnedfilename, 
                        const char* pcExceptionfilename,
                        const char* pcfilename);

private:
    bool Loadbinfile(const char* filename);
    bool LoadExceptionfile(const char* filename);
    int GetIndex(char cChromosome, int nOffset);
    void SearchRecords(size_t index,char cChrnum,size_t qstart, size_t qstop);
    void GetExceptionRecord(char cChar, size_t qstart, size_t qstop);
    
    stAnnotateRecs m_Records;    
    stBinnedIdx m_aBinnedIdx[24];
    stExceptionIdx m_aExceptionIdx[24];
    ifstream m_snpfile;
    char *m_pcAnnotate;
    int m_nfdAnnotate;
    size_t m_file_size;    
};

#endif	/* _SXANNOTATE_H */

