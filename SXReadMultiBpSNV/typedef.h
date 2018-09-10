#include <iostream>
using namespace std;

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
} stSNV;
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


typedef enum eSNVType{eSNP='S',eINS='I',eDEL='D', eALL='A'};

////////////////////////////////////////////////////////////////////////////////

typedef struct _stSNVFile{
   int fd;
   unsigned long ulFileSize;
   stSNV *pstSNV;

   _stSNVFile()
   {
        fd=0; ulFileSize=0; pstSNV=NULL;
   }
} stSNVFile;

////////////////////////////////////////////////////////////////////////////////

const unsigned int SNVREC_SIZE = 1000000;
const unsigned long INIT_MERGED_LIST_SIZE = 500000000;
const unsigned long EXT_MERGED_LIST_SIZE = 100000000;
   

typedef struct _stSNVRec{
    //unsigned char acOffset[4];
    unsigned int unOffset;   
    char cRefBase;
    char cVarBase;
    unsigned char ucQryPos;
    unsigned char ucQltyScore;
    unsigned short usPEnd;
    char cStrand;
}stSNVRec;


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


typedef struct _stSNVClusterRec{
    unsigned int unOffset;
    unsigned int unCnt;
    unsigned int unFwdCnt;
    unsigned int unRvsCnt;
    unsigned int unPEndCnt;
    unsigned int ulTotalQS;
    unsigned int ulTotalQP;     
    unsigned char *pucMisc;
    unsigned char ucMaxQS;  
    char cChromosome;  
    char cRefBase;
    char cVarBase;
    //bool bKeep;

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

 
    void Add(unsigned char **pucList, unsigned char ucChar)
    {         
    
         if (*pucList == NULL){
             *pucList = (unsigned char*)malloc(2); 
             (*pucList)[0]=ucChar; (*pucList)[1]='\0';
         }  
            
         else{            
             int nLen=strlen((const char*)(*pucList))+2;  
             (*pucList) = (unsigned char*)realloc((*pucList),nLen);
             (*pucList)[nLen-2]= ucChar; (*pucList)[nLen-1]='\0';
         }
    }


   void AddMisc(unsigned char ucSide,unsigned char ucStrand,unsigned char ucQryPos)
   {
      if (!pucMisc){
        pucMisc = (unsigned char*)malloc(6);
        sprintf((char*)&pucMisc[0],"%c%c%u%c",ucSide,ucStrand, ucQryPos,'\0');
      } 
      else{
        int nLen=strlen((const char*)(pucMisc))+7; 
        pucMisc = (unsigned char*)realloc(pucMisc,nLen);
        sprintf((char*)&pucMisc[nLen-7],",%c%c%u%c",ucSide,ucStrand, ucQryPos,'\0');
      }       

      
   } 
   
}stSNVClusterRec;

