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

typedef enum eSNVType{eSNP='S',eINS='I',eDEL='D', eALL='A'};

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
    char cChromosome;
    unsigned char acOffset[4];
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand; //0-rvs, 1-fwd
    unsigned long ulLnkFileOffset;//char acLnkFileOffset[8];
} stSNV;

#pragma pack()

typedef struct _compare{
    bool operator()(const stSNV *st1,const stSNV *st2)
    {
        if (st1->ucQltyScore != st2->ucQltyScore)
            return st1->ucQltyScore < st2->ucQltyScore;

        return st1->ucQryPos < st2->ucQryPos;
    }
};


typedef std::map<stSNV*, unsigned long, _compare> SNV_STATS_LIST;


typedef struct _compareChromosome{
    bool operator()(const stSNV *st1, const stSNV *st2)
    {
        return abs(st1->cChromosome) < abs(st2->cChromosome);
    }
};

typedef std::map<stSNV*, unsigned long,_compareChromosome> CHRO_FREQ_LIST;


typedef struct _stSNVFreq
{
    char cChromosome;
    unsigned int unOffset;
    char cStrand;
    char cRefBase;
    char cVarBase;

    _stSNVFreq(stSNV *pSNV)
    {
        cChromosome = abs(pSNV->cChromosome);
        unOffset = ComputeOffset(pSNV->acOffset);
        cStrand = (pSNV->cChromosome > 0)?'+':'-';
        cRefBase = pSNV->cRefBase;
        cVarBase = pSNV->cVarBase;
    }
} stSNVFreq;


typedef struct _compareFreq{
    bool operator()(const stSNVFreq *st1, const stSNVFreq *st2)
    {
        if (st1->cChromosome != st2->cChromosome)
            return st1->cChromosome < st2->cChromosome;

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
    char cChromosome;
    unsigned int unOffset;

    _stChromOffset(stSNV *pSNV)
    {
        cChromosome = abs(pSNV->cChromosome);
        unOffset = ComputeOffset(pSNV->acOffset);
    }

    _stChromOffset(char chromo, unsigned int Offset)
    {
        cChromosome = abs(chromo);
        unOffset = Offset;
    }
} stChromOffset;


typedef struct _compareChromOffset{
    bool operator()(const stChromOffset *st1, const stChromOffset *st2)
    {
        if (st1->cChromosome != st2->cChromosome)
            return st1->cChromosome < st2->cChromosome;

        return st1->unOffset < st2->unOffset;
    }
};

typedef std::map<stChromOffset*, unsigned long, _compareChromOffset> CHRO_OFFSET_FREQ_LIST;


typedef struct _stChroStrand
{
    char cChromosome;
    char cStrand;

    _stChroStrand(stSNV *pSNV)
    {
        cChromosome = abs(pSNV->cChromosome);
        cStrand = (pSNV->cChromosome > 0)?'+':'-';
    }
}stChroStrand;


typedef struct _compareChroStrand{
    bool operator()(const stChroStrand *st1, const stChroStrand *st2)
    {
        if (st1->cChromosome != st2->cChromosome)
            return st1->cChromosome < st2->cChromosome;

        return st1->cStrand < st2->cStrand;
    }
};

typedef std::map<stChroStrand*, unsigned long, _compareChroStrand> CHRO_STRAND_LIST;


typedef struct _compareBase{
    bool operator()(const char c1, const char c2)
    {
        return c1 < c2;
    }
};

////////////////////////////////////////////////////////////////////////////////

#pragma pack(1)

typedef struct _stSNVRec{
    unsigned char acOffset[4];
    unsigned int unOffset;
    //char cChromosome;
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand;
    unsigned long ulLnkFileOffset;
}stSNVRec;


typedef struct _stSNVClusterRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int unFwdCnt;
    unsigned int unRvsCnt; 
    unsigned int unPEndCnt;
    unsigned long ulTotalQS;
    unsigned long ulTotalQP;    
    unsigned char *pucMisc; 
    unsigned char ucMaxQS;     
    char cChromosome;
    char cRefBase;
    char cVarBase;
    bool bKeep;

    _stSNVClusterRec()
    {
        InitRec();    
    }

    ~_stSNVClusterRec()
    {
        ClrRec();
    }

    void InitRec()
    {
        //bKeep = true;
        unCnt=unFwdCnt=unRvsCnt=unPEndCnt=0;
        ulTotalQP=ulTotalQS=0; ucMaxQS=0;
        pucMisc = NULL;
    }

    void ClrRec()
    {
       if (pucMisc) free(pucMisc);
    }


    void AddMisc(const unsigned char &ucSide, const unsigned char &ucStrand,
                const unsigned char &ucQryPos, const char &cPEType,  
                const char &cAmbiguity, const unsigned long &ulLnkFileOffset)
    {
       if (!pucMisc){
           pucMisc = (unsigned char*)malloc(28);
           sprintf((char*)&pucMisc[0],"%c%c%u%c%c%lu%c",ucSide,ucStrand,ucQryPos,cPEType,cAmbiguity,ulLnkFileOffset,'\0');
       }
       else{
           int nLen=strlen((const char*)(pucMisc))+29;
           pucMisc = (unsigned char*)realloc(pucMisc,nLen);
           sprintf((char*)&pucMisc[nLen-29],",%c%c%u%c%c%lu%c",ucSide,ucStrand, ucQryPos,cPEType,cAmbiguity,ulLnkFileOffset,'\0');
       }
    }

}stSNVClusterRec;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct _stMBSNVRec{
    unsigned int unOffset;
    unsigned char ucBpsCnt;     
    unsigned char *pucMisc; //points to Refbase, VarBase,QScore
    unsigned char ucQryPos;
    unsigned short usPEnd;
    char cStrand;  //Rightmost bit-Strand, SNVType(2bits):0-SNP,1-INS,2-DEL
    unsigned long ulLnkFileOffset;
}stMBSNVRec;


typedef struct _stMBSNVClusterRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int unFwdCnt;
    unsigned int unRvsCnt;
    unsigned int unPEndCnt;
    //unsigned long ulTotalQS;
    //unsigned long ulTotalQP;
    unsigned char *pucMisc;
    //unsigned char ucMaxQS;
    unsigned char ucChromosome;//Rightmost 5 bits:Chro, follow by:SNVType,bKeep
    unsigned char ucBpsCnt;
    //char cVarBase;
    //bool bKeep;

    _stMBSNVClusterRec()
    {
        InitRec();
    }

    ~_stMBSNVClusterRec()
    {
        ClrRec();
    }

    void InitRec()
    {
        //bKeep = true;
        unCnt=unFwdCnt=unRvsCnt=unPEndCnt=0;
        //ulTotalQP=ulTotalQS=0; ucMaxQS=0;
        pucMisc = NULL;
    }

    void ClrRec()
    {
       if (pucMisc) free(pucMisc);
    }

    void AddMisc(const unsigned char &ucSide, const unsigned char &ucStrand,
                const unsigned char &ucQryPos, const char &cPEType,
                const char &cAmbiguity, const unsigned long &ulLnkFileOffset)
    {
         int nLen=strlen((const char*)(pucMisc))+28;
         pucMisc = (unsigned char*)realloc(pucMisc,nLen);
         sprintf((char*)&pucMisc[nLen-28],"%c%c%u%c%c%lu%c",ucSide,ucStrand,ucQryPos,cPEType,cAmbiguity,ulLnkFileOffset,'\0');
    }

    void InitMisc(const unsigned char *pcVal, int nBpsCnt, const char &cStrand)
    {
         if (pucMisc) {free(pucMisc);}
         int nLen, nNumOfQScore, nTotalAllele = nBpsCnt*2;

         nNumOfQScore = ((cStrand & 0x06)== 0x04/*DEL*/)?1:nBpsCnt;
         nLen = nTotalAllele+(nNumOfQScore*3); 
         pucMisc = (unsigned char*)malloc(nLen+1); pucMisc[nLen]='\0';                    
         memcpy((char*)&pucMisc[0],pcVal,nTotalAllele);

         for (unsigned int i=nTotalAllele,j=i; i<strlen((const char*)pcVal); i++,j+=3) 
		 sprintf((char*)&pucMisc[j],"%03u",pcVal[i]);           

         pucMisc[nLen]='\0';
    }

    void AddQScore2Misc(const unsigned char *pcVal, int nBpsCnt, const char &cStrand)
    {
         int nOrgLen, nLen, nNumOfQScore, nTotalAllele = nBpsCnt*2;

         nNumOfQScore = ((cStrand & 0x06)== 0x04/*DEL*/)?1:nBpsCnt;  
         nOrgLen = strlen((const char*)pucMisc);
         nLen = nOrgLen + (nNumOfQScore*3)+1; 

         pucMisc = (unsigned char*)realloc(pucMisc,nLen+1); pucMisc[nLen]=0;
         pucMisc[nOrgLen] = ',';

         for (unsigned int i=nTotalAllele,j=nOrgLen+1; i<strlen((const char*)pcVal); i++, j+=3)
             sprintf((char*)&pucMisc[j],"%03u",pcVal[i]);            
    }
   
}stMBSNVClusterRec;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

typedef struct _stSNPRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int ulTotalQS;
    char cChromosome;
    char cRefBase;
    char cVarBase;
}stSNPRec;

typedef struct _stINSRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int ulTotalQS;
    char cChromosome;
    char cVarBase;
}stINSRec;


typedef struct _stDELRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int ulTotalQS;
    char cChromosome;
    char cRefBase;
}stDELRec;

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
    //unsigned long ulTotalQS;
    unsigned long ulTotalReadDens;

    _SNVTblStats()
    {
        //ulTotalQS=
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

typedef struct _CompareGeneKey
{
    bool operator()(const char* pc1, const char *pc2)
    {
        return strcmp(pc1,pc2)<0;
    }
};

typedef std::map<const char*,const char*,_CompareGeneKey> GENE_NAME_LIST;

////////////////////////////////////////////////////////////////////////////////

typedef struct _SNVAnalysis{
    unsigned long ulTotalSNP;
    unsigned long ulTotalINS;
    unsigned long ulTotalDEL;
    unsigned long ulTotalSNP_dbsnp;
    unsigned long ulTotalINS_dbsnp;
    unsigned long ulTotalDEL_dbsnp;
    unsigned long ulTotalSNP_gene;
    unsigned long ulTotalSNP_exon;
    unsigned long ulTotalINS_exon;
    unsigned long ulTotalDEL_exon;
    unsigned long ulTotalSNP_Transitions;
    unsigned long ulTotalSNP_Transversions;

    _SNVAnalysis()
    {
        ulTotalSNP_gene=ulTotalSNP=ulTotalINS=ulTotalDEL=ulTotalSNP_dbsnp=
        ulTotalINS_dbsnp=ulTotalDEL_dbsnp=ulTotalSNP_exon=
        ulTotalINS_exon=ulTotalDEL_exon=ulTotalSNP_Transitions=
        ulTotalSNP_Transversions=0;
    }
}stSNVAnalysis;


//typedef enum eSNVType{eSNP='S',eINS='I',eDEL='D', eALL='A'};

////////////////////////////////////////////////////////////////////////////////

typedef struct _compareQryPos{
    bool operator()(unsigned char p1, unsigned char p2)
    {
        return p1 < p2;
    }
};

typedef std::map<unsigned char, unsigned long,_compareQryPos> QRY_POS_LIST;


typedef struct _compareQlyScore{
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

typedef struct _stSNVTrace
{
    char cChromosome;
    unsigned int unOffset;
    char cRefBase;
    char cVarBase;

    _stSNVTrace(stSNV *pSNV)
    {
        cChromosome = abs(pSNV->cChromosome);
        unOffset = ComputeOffset(pSNV->acOffset);
        cRefBase = pSNV->cRefBase;
        cVarBase = pSNV->cVarBase;
    }
}stSNVTrace;


typedef struct _compareSNVTrace{
    bool operator()(const stSNVTrace *st1, const stSNVTrace *st2)
    {
        if (st1->cChromosome != st2->cChromosome)
            return st1->cChromosome < st2->cChromosome;

        if (st1->unOffset != st2->unOffset)
            return st1->unOffset < st2->unOffset;

        if (st1->cRefBase != st2->cRefBase)
            return st1->cRefBase < st2->cRefBase;

        return st1->cVarBase < st2->cVarBase;
    }
};

typedef std::map<stSNVTrace*,stSNVTraceData*,_compareSNVTrace>SNV_TRACE_LIST;


typedef struct _CompareGeneKeyword
{
    bool operator()(const char* pc1, const char *pc2)
    {
        return strcmp(pc1,pc2)<0;
    }
};

typedef std::map<const char*,const char*,_CompareGeneKeyword> GENE_KEYWORD_LIST;
typedef std::vector <char*> VEC_LIST;

#endif	/* _TYPEDEF_H */

