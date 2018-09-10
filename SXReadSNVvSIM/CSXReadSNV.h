/* 
 * File:   CSXReadSNV.h
 * Author: Hwah Shih Yiew
 *
 * Created on November 17, 2009, 12:40 PM
 */


#ifndef _CSXREADSNV_H
#define	_CSXREADSNV_H

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include "typedef.h"
#include "SXAnnotate.h"


#define SNVREC_SIZE 1000000

typedef struct _stSNVRecList
{
    stSNVRec *pSNVRec;
    int nTotal;
    int nUsed;
    pthread_t ThreadID;

    _stSNVRecList()
   { 
       pSNVRec = (stSNVRec*)malloc(sizeof(stSNVRec)*SNVREC_SIZE);
       nUsed = 0;
       nTotal = SNVREC_SIZE; 
   }

   ~ _stSNVRecList()
   {
       if (pSNVRec) free(pSNVRec);
   }

}stSNVRecList;



class CSXReadSNV {
public:
    CSXReadSNV();
    virtual ~CSXReadSNV();

    void run();
    void setQltyScore(int nQltyScore1,int nQltyScore2,const char* pcFile);
    void setQryPos(int nQryPos1,int nQryPos2,const char* pcFile);
    void setFilterParams(int nQryPos1,int nQryPos2,int nQltyScore1,
                         int nQltyScore2,char cSNVType,const char* pcFile);
    
    void setStateRpt(const char* pcFile);
    void setSNPStateRpt(const char* pcFile);
    void setINSStateRpt(const char* pcFile);
    void setDELStateRpt(const char* pcFile);
    void setAvgQSLstInputs(const char *pcdbSnp, const char *pcdbIndel,
                           const char *pcDensFilePath, unsigned int unMinSupportingReads,
                           unsigned int unMaxSupportingReads,float fMinReadStrength,
                           int nMinQScore, const char *pcSampleID, const char *pcGenomeVer,
                           int nQScoreOfs);
    void setSNVLstInputs(const char *pcBinnedFile,const char *pcExceptionFile,
                         const char *pcAnnotateFile, unsigned int unMinSupportingReads,
                         unsigned int unMaxSupportingReads, char cSNVType,
                         const char *pcSampleID, const char *pcGenomeVer);
    
    void setSNVTraceInputs(unsigned int unSupportingReads,char cSNV);
    void setDensityInput(const char* pcFile, char cSNVType, const char* pcSampleID);
    void OutputREADDensity1K();
    

    bool m_bOut2Screen;
    bool m_bFilter;
    bool m_bFilterx;
    bool m_bSNVTrace;

    const char *m_pcInFile;
    const char *m_pcOutFile;
    const char *m_pcNewFile;
    
    const char *m_pcFreqFile;
    const char *m_pcChroFreqFile;
    const char *m_pcChroStrandFile;
    const char *m_pcChromOffsetFile;
    const char *m_pcStatsFile;
    const char *m_pcSNPStatsFile;
    const char *m_pcINSStatsFile;
    const char *m_pcDELStatsFile;    
    const char *m_pcDensFilePath;
    const char *m_pcSampleID;
    const char *m_pcGenomeVer;

    unsigned int m_unMinSupportingReads;

private:
    void GetInitParam();
    void FilterSNVList2FileByQP();
    void FilterSNVList2FileByQS();
    void FilterSNVList();
    void FilterSNVList_SNP();
    void FilterSNVList_INS();
    void FilterSNVList_DEL();
    void FilterSNVListx();
    void FilterSNVListx_SNP();
    void FilterSNVListx_INS();
    void FilterSNVListx_DEL();
    void ClrMemoryMap();
    void OutputInVal2Scrn();
    void OutputInVal2File();
    void DisplaySNVStats();
    void OutputSNVStats(const char c, const char* pcFile);
    void OutputSNVFreq();
    void OutputChroFreq();
    void OutputChroStrand();
    void OutputChromOffset();
    void OutputDensity();
    void OutputSNVList();
    void PrintSNPList(FILE *pfSrc);
    void PrintINSList(FILE* pfSrc);
    void PrintDELList(FILE* pfSrc);
    void PrintStatsRpt(FILE *pf, stSNVTblStats *pstTblStats,char *pcType);
    void PrintSNVAnalysis();
    void OutputType(FILE *pf, char cRefBase, char cVarBase);
    bool LoadGeneName(FILE *pf);
    void ClrGeneNameMap();
    void TraceSNVSupportingReads();
    void OutputAvgQScoreTable();
    void PrintAvgQScoreList();
    void PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType);
    void FilterSNVRec(stSNVClusterRec &Rec1, stSNVClusterRec &Rec2);
    void FilterSNVRec(stSNVClusterRec &Rec1, stSNVClusterRec &Rec2,
                      stSNVClusterRec &Rec3);
    void PrintBaseRepeat(FILE *pf,stSNVClusterRec &Rec);
    void LogInvalidRec(FILE *pf, stSNVRec &Rec, int nChromosome);

    inline void GenerateSNVTbl();
    inline void FilterSNVTbl();

   
    inline double CalcAvgQS(unsigned long ulTotalQS, unsigned long ulTotalCnt);
    inline void OutputChrom(FILE *pf, unsigned int unChrom);
    inline unsigned int GetReadDensity(FILE *pf, unsigned int uioffset);
    inline unsigned int GetReadDensityEx(FILE *pf, unsigned int uioffset);
    inline double Percentage(unsigned long ulVal, unsigned long ulDenominator);
    inline double Division(unsigned long ulTotal, unsigned long ulDenominator);
    inline void PrintRptDateTime(FILE *pf);
    inline void AllocMergeListSize();
    inline void ExtractAnnotation(char *pRecs, stAnnotate &aAnnotate);
    inline void PrintVECListContents(VEC_LIST &vec, FILE *pf);

        
    stSNV *m_pSNV;
    stSNVRec *m_pSNVRec;
    stSNPRec *m_pSNP;
    stSNVClusterRec *m_pSNVClustered;
    SNV_STATS_LIST m_SNVStats_List;   
    SNV_FREQ_LIST m_SNVFreq_List;
    CHRO_FREQ_LIST m_ChroFreq_List;    
    CHRO_STRAND_LIST m_ChroStrand_List;
    stSNVTblResult m_stSNVTblResult;
    stSNVAnalysis m_stSNVAnalysis;
    SXAnnotate *m_pSXAnnotate;


    FILE **m_apfRD;

    size_t m_file_size;

    eSNVType m_eSNV;

    unsigned int m_unMaxSupportingReads;
    unsigned int m_unMaxReadDensity;
    unsigned int m_unMinBinDev;
    unsigned long m_ulTotalRecs;
    unsigned long INIT_MERGED_LIST_SIZE;
    unsigned long EXT_MERGED_LIST_SIZE;

    int m_nFQryPos1;
    int m_nFQryPos2;
    int m_nFQltyScore1;
    int m_nFQltyScore2;
    int m_nQScoreOfs;

    bool m_bSNVLst;
    bool m_bAvgQScoreLst;
    bool m_bMemMapped;    
    bool m_bFQryPos;
    bool m_bFQltyScore;
    bool m_bStateRpt;
    bool m_bSNPStateRpt;
    bool m_bINSStateRpt;
    bool m_bDELStateRpt;
    bool m_fwdflag;
    bool m_revflag;
    bool m_bDensity;

    const char *m_pcSNVTbl;
    const char *m_pcBinnedFile;
    const char *m_pcExceptionFile;
    const char *m_pcAnnotateFile;
    const char *m_pcdbSnp;
    const char *m_pcdbIndel;
    
    const char *m_pcCNV;
    
    //float m_fAvgQScore;
    float m_fMinReadStrength;
    int m_nMinQScore;

    size_t merge_list_alloc, merge_list_used;
    unsigned int *m_punGINums;
};

#endif	/* _CSXREADSNV_H */

