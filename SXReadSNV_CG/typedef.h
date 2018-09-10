/* 
 * File:   typedef.h
 * Author: Hwah Shih Yiew
 *
 * Created on November 19, 2009, 4:41 PM
 */

#ifndef _TYPEDEF_H
#define	_TYPEDEF_H

#include <map>
#include <vector>


inline unsigned int ComputeOffset(unsigned char acVal[4])
{
    unsigned int unOffset= acVal[3];unOffset<<=8;
    unOffset|= acVal[2];unOffset<<=8;
    unOffset|= acVal[1];unOffset<<=8;
    unOffset|= acVal[0];

    return unOffset;
}

#pragma pack(1)
typedef struct _stSNV
{
    int nChromosome;
    unsigned char acOffset[4];
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand; //0-rvs, 1-fwd
} stSNV;

#pragma pack()

struct _compare{
    bool operator()(const stSNV *st1,const stSNV *st2)
    {
        if (st1->ucQltyScore != st2->ucQltyScore)
            return st1->ucQltyScore < st2->ucQltyScore;

        return st1->ucQryPos < st2->ucQryPos;
    }
};


typedef std::map<stSNV*, unsigned long, _compare> SNV_STATS_LIST;


struct _compareChromosome{
    bool operator()(const stSNV *st1, const stSNV *st2)
    {
        return st1->nChromosome < st2->nChromosome;
    }
};

typedef std::map<stSNV*, unsigned long,_compareChromosome> CHRO_FREQ_LIST;


typedef struct _stSNVFreq
{
    int nChromosome;
    unsigned int unOffset;
    char cStrand;
    char cRefBase;
    char cVarBase;

    _stSNVFreq(stSNV *pSNV)
    {
        nChromosome = pSNV->nChromosome;
        unOffset = ComputeOffset(pSNV->acOffset);
        cStrand = (pSNV->cStrand > 0)?'+':'-';
        cRefBase = pSNV->cRefBase;
        cVarBase = pSNV->cVarBase;
    }
} stSNVFreq;


struct _compareFreq{
    bool operator()(const stSNVFreq *st1, const stSNVFreq *st2)
    {
        if (st1->nChromosome != st2->nChromosome)
            return st1->nChromosome < st2->nChromosome;

        if (st1->unOffset != st2->unOffset)
            return st1->unOffset < st2->unOffset;

        if (st1->cStrand != st2->cStrand)
            return st1->cStrand > st2->cStrand;

        if (st1->cRefBase != st2->cRefBase)
            return st1->cRefBase < st2->cRefBase;

        return st1->cVarBase < st2->cVarBase;
    }
};

typedef std::map<stSNVFreq*, unsigned long, _compareFreq> SNV_FREQ_LIST;


typedef struct _stChromOffset
{
    int nChromosome;
    unsigned int unOffset;

    _stChromOffset(stSNV *pSNV)
    {
        nChromosome = pSNV->nChromosome;
        unOffset = ComputeOffset(pSNV->acOffset);
    }

    _stChromOffset(int nChromo, unsigned int Offset)
    {
        nChromosome = nChromo;
        unOffset = Offset;
    }
} stChromOffset;


struct _compareChromOffset{
    bool operator()(const stChromOffset *st1, const stChromOffset *st2)
    {
        if (st1->nChromosome != st2->nChromosome)
            return st1->nChromosome < st2->nChromosome;

        return st1->unOffset < st2->unOffset;
    }
};

typedef std::map<stChromOffset*, unsigned long, _compareChromOffset> CHRO_OFFSET_FREQ_LIST;


typedef struct _stChroStrand
{
    int nChromosome;
    char cStrand;

    _stChroStrand(stSNV *pSNV)
    {
        nChromosome = pSNV->nChromosome;
        cStrand = (pSNV->cStrand > 0)?'+':'-';
    }
}stChroStrand;


struct _compareChroStrand{
    bool operator()(const stChroStrand *st1, const stChroStrand *st2)
    {
        if (st1->nChromosome != st2->nChromosome)
            return st1->nChromosome < st2->nChromosome;

        return st1->cStrand < st2->cStrand;
    }
};

typedef std::map<stChroStrand*, unsigned long, _compareChroStrand> CHRO_STRAND_LIST;

////////////////////////////////////////////////////////////////////////////////

#pragma pack(1)

typedef struct _stSNVRec{
    unsigned char acOffset[4];
    unsigned int unOffset;
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand;
}stSNVRec;


typedef struct _stSNVClusterRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int unFwdCnt;
    unsigned int unRvsCnt; 
    unsigned int unPEndCnt;
    unsigned long ulTotalQS;
    unsigned long ulTotalQP;    
    unsigned char ucMaxQS;
    int nChromosome;
    char cRefBase;
    char cVarBase;
    bool bKeep;

    _stSNVClusterRec()
    {
        bKeep = true;
        unCnt=unFwdCnt=unRvsCnt=0;
        ulTotalQP=0; ulTotalQS=0; ucMaxQS=0; unPEndCnt=0;
    }
}stSNVClusterRec;

#pragma pack()
////////////////////////////////////////////////////////////////////////////////

typedef struct _SNVTblResult{
    unsigned int unACGT[4];
    unsigned int unDEL;

    unsigned long ulTotalQS[4];
    unsigned long ulDelTotalQS;


    _SNVTblResult()
    {
        for (short i=0; i<4;i++)
        {
            unACGT[i]=ulTotalQS[i]=0;
        }

        unDEL=ulDelTotalQS=0;
    }

}stSNVTblResult;


typedef struct _SNVTblStats{
    unsigned int unCnt;
    unsigned int undbSNP;
    unsigned int unGene;
    unsigned int unExon;
    unsigned long ulTotalReadDens;

    _SNVTblStats()
    {
        unCnt=undbSNP=unGene=unExon=ulTotalReadDens=0;

    }
}stSNVTblStats;


typedef struct _stAnnotate
{
    std::vector<char*> GeneName_List;
    std::vector <char*> GeneDesc_List;
    std::vector <char*> GeneKW_List;
    std::vector<char*> ExonID_List;
    std::vector <char*> CNV_List;
    std::vector <char*> miRNA_List;
    std::vector <char*> Promoter_List;
    std::vector <char*> UTR_List;

    ~ _stAnnotate()
    {
        ClrList();
    }

    void ClrList()
    {
       GeneName_List.clear(); GeneDesc_List.clear();
       GeneKW_List.clear(); ExonID_List.clear();
       CNV_List.clear(); miRNA_List.clear(); Promoter_List.clear();
       UTR_List.clear();  
    }

} stAnnotate;

////////////////////////////////////////////////////////////////////////////////

typedef struct _SNVAnalysis{
    unsigned long ulTotalSNP;
    //unsigned long ulTotalINS;
    //unsigned long ulTotalDEL;
    unsigned long ulTotalSNP_dbsnp;
    //unsigned long ulTotalINS_dbsnp;
    //unsigned long ulTotalDEL_dbsnp;
    unsigned long ulTotalSNP_gene;
    unsigned long ulTotalSNP_exon;
    //unsigned long ulTotalINS_exon;
    //unsigned long ulTotalDEL_exon;
    unsigned long ulTotalSNPHet;
    unsigned long ulTotalSNPHom;
    unsigned long ulTotalSNPGreyArea; 
    unsigned long ulTotalSNP_Transitions;
    unsigned long ulTotalSNP_Transversions;

    _SNVAnalysis()
    {
        /*ulTotalSNP_gene=ulTotalSNP=ulTotalINS=ulTotalDEL=ulTotalSNP_dbsnp=
        ulTotalINS_dbsnp=ulTotalDEL_dbsnp=ulTotalSNP_exon=
        ulTotalINS_exon=ulTotalDEL_exon=ulTotalHet=ulTotalHom=
        ulTotalGreyArea=ulTotalSNP_Transitions=ulTotalSNP_Transversions=0;
        */

        Init(); 
    }

   void Init()
   {
        ulTotalSNP_gene=ulTotalSNP=ulTotalSNP_dbsnp=ulTotalSNP_exon=ulTotalSNPHet=
        ulTotalSNPHom=ulTotalSNPGreyArea=ulTotalSNP_Transitions=ulTotalSNP_Transversions=0;
   }  
}stSNVAnalysis;


enum eSNVType{eSNP='S',eINS='I',eDEL='D', eALL='A'};

////////////////////////////////////////////////////////////////////////////////

struct _compareQryPos{
    bool operator()(unsigned char p1, unsigned char p2)
    {
        return p1 < p2;
    }
};

typedef std::map<unsigned char, unsigned long,_compareQryPos> QRY_POS_LIST;


struct _compareQlyScore{
    bool operator()(unsigned char p1, unsigned char p2)
    {
        return p1 < p2;
    }
};

typedef std::map<char, unsigned long, _compareQlyScore> QLY_SCORE_LIST;

typedef struct _stSNVTraceData
{
    unsigned int unCnt;
    QRY_POS_LIST QryPosList;
    QLY_SCORE_LIST QlyScoreList;
}stSNVTraceData;

typedef std::vector <char*> VEC_LIST;

#endif	/* _TYPEDEF_H */

