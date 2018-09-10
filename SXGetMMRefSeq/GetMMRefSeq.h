/* 
 * File:   GetMMRefSeq.h
 * Author: Hwah Shih Yiew
 *
 * Created on April 13, 2010, 12:14 PM
 */

#ifndef _GETHREFSEQ_H
#define	_GETHREFSEQ_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>

class GetMMRefSeq {
public:
    GetMMRefSeq();
    virtual ~GetMMRefSeq();

    void run_SV(const char *pcInput, const char *pcRefFilePath,
                const char *pcOutput, int nSize);

    void run_SnpIndel(const char *pcInput, const char *pcRefFilePath,
                      const char *pcOutput, int nSize);

private:

    FILE *m_pfIn, *m_pfOut;
    FILE **m_apfMMRef;
    char **m_apcMMRef;
    int m_anMMRefSize[24];

    void LoadMMRefSeq();
    bool OpenFiles(const char *pcInput, const char *pcRefFilePath,
                   const char *pcOutput);
    inline void CloseRefFiles();
    void CloseAllFiles();
};

#endif	/* _GETMMREFSEQ_H */

