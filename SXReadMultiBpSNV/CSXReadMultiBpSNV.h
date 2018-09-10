#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "typedef.h"



class CSXReadMultiBpSNV{
public:
    CSXReadMultiBpSNV();
    virtual ~CSXReadMultiBpSNV();

    void PrintSNVFileContents(char *pcOutput);
    void setMbpSNVListInput(const char *pcdbSnp, const char *pcdbIndel,
                            const char *pcDensFilePath,
                            unsigned int unMinSupportingReads,
                            unsigned int unMaxSupportingReads,
                            float fMinReadStrength,int nMinQScore,
                            const char *pcSampleID,const char *pcGenomeVer); 

    void Cluster(); 
    const char *m_pcInFile;

private:
    void ClrMemoryMap();
    bool MapSNVFile();
    void GenerateSNVTbl();
    void AllocMergeListSize();
    void PrintMultiBpSNVList(); 
    //void Cluster(); 
    void PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType); 
    void CSXReadMultiBpSNV::SetPEndnStrandInfo(stSNVClusterRec&, short unsigned int&, char&);
    void SetMiscInfo(stSNVClusterRec &SNVClusterRec, unsigned short &usPEnd, char &cStrand, 
                     unsigned short &usQryPos);

    inline double CalcAvgQS(unsigned long ulTotalQS, unsigned long ulTotalCnt);
    inline void PrintRptDateTime(FILE *pf);
    inline void PrintRemarksParam(FILE *pfParam); 
    inline void OutputChrom(FILE *pf, unsigned int unChrom); 
    inline unsigned int GetReadDensity(FILE *pf, unsigned int uioffset);
    inline unsigned int GetReadDensityEx(FILE *pf, unsigned int uioffset);  
    inline double Percentage(unsigned long ulVal, unsigned long ulDenominator);
    inline double Division(unsigned long ulTotal, unsigned long ulDenominator);

    FILE **m_apfRD; 
    stSNV *m_pSNV;
    stSNVFile m_stSNVFile;
    stSNVClusterRec *m_pSNVClustered;
    stSNVAnalysis m_stSNVAnalysis; 

    const char *m_pcdbSnp; 
    const char *m_pcdbIndel;
    const char *m_pcDensFilePath;

    unsigned int m_unMinSupportingReads;
    unsigned int m_unMaxSupportingReads;
    float m_fMinReadStrength; 
    int m_nMinQScore;
    const char *m_pcSampleID;
    const char *m_pcGenomeVer;
    bool m_bMultiBpSNVList;

    int m_nQScoreOfs;
    unsigned long merge_list_alloc;
    unsigned long merge_list_used;
};

