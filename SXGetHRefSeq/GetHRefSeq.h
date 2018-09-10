/* 
 * File:   GetHRefSeq.h
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

class GetHRefSeq {
public:
    GetHRefSeq();
    virtual ~GetHRefSeq();

    void run_SV(const char *pcInput, const char *pcRefFilePath,
                const char *pcOutput, int nSize);

    void run_SnpIndel(const char *pcInput, const char *pcRefFilePath,
                      const char *pcOutput, int nSize);

private:

    FILE *m_pfIn, *m_pfOut;
    FILE **m_apfHsRef;
    char **m_apcHsRef;
    int m_anHsRefSize[24];

    void LoadHRefSeq();
    bool OpenFiles(const char *pcInput, const char *pcRefFilePath,
                   const char *pcOutput);
    inline void CloseRefFiles();
    void CloseAllFiles();
};

#endif	/* _GETHREFSEQ_H */

