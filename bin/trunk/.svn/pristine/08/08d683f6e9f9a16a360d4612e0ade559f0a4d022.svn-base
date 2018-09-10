#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const int REC_SIZE = 1000000;

typedef struct _stRec
{   
    int nOffset;
    char acAllele[4];
    float fCfdScore;

    _stRec()
    {
        nOffset=0; fCfdScore=0.0;
    }
}stRec;


typedef struct _stSNVList
{
    FILE *pf;    
    int anTotal[24]; 
    int anUsed[24];
    stRec *apstRec[24];

    _stSNVList()
    {
        pf = NULL; 
        for (int i=0; i<24; i++){  
            apstRec[i] = (stRec*)malloc(sizeof(stRec)*REC_SIZE);
            anUsed[i]=0; anTotal[i]=REC_SIZE;
        }
    }

    ~_stSNVList()
    {
        if (pf) fclose(pf);
        for (int i=0; i<24; i++){
            if (apstRec[i]) free(apstRec[i]);
        }  
    }
}stSNVList;

const int SEQSNV_SIZE = 120;

typedef struct _stSeqSNV
{
    char *pcData;
    int nTotal; 
    int nUsed;

    _stSeqSNV()
    {
        pcData = (char*)malloc(SEQSNV_SIZE);
        nUsed=0; nTotal=SEQSNV_SIZE; 
    }  

    ~_stSeqSNV()
    {
       if (pcData) free(pcData);  
    } 

   void Add(char *pcIn, int nLen)
    {               
        if (nTotal < (nUsed+nLen)) {
            nTotal += (SEQSNV_SIZE+nLen);
            pcData = (char*)realloc(pcData,nTotal);
        } 

        memcpy(pcData+nUsed,pcIn,nLen); nUsed += nLen;  
    } 

   void Add(char Chr)
   {
      if (nUsed == nTotal){
         nTotal += SEQSNV_SIZE;
         pcData = (char*)realloc(pcData,nTotal);  
      }      
      pcData[nUsed]=Chr; nUsed++;
   }
}stSeqSNV;



typedef enum eSNVType{eSNP=0,eINS,eDEL};

FILE *g_pfSeq=NULL,*g_pfOut=NULL;
int g_nSize;
stSNVList g_stSNVList[3];
char g_cType;

