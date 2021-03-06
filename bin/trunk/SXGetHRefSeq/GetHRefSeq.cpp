/* 
 * File:   GetHRefSeq.cpp
 * Author: Manish
 * 
 * Created on April 13, 2010, 12:14 PM
 */

#include "GetHRefSeq.h"

#define _TESTING 0

GetHRefSeq::GetHRefSeq() {

    m_pfIn=NULL;m_pfOut=NULL;

    m_apfHsRef = new FILE*[24]; m_apcHsRef = new char*[24];
    
    for (short i=0;i<24;i++){
        m_apfHsRef[i]=NULL; m_apcHsRef[i]=NULL;
    }

}


GetHRefSeq::~GetHRefSeq() {
    if (m_pfIn) fclose(m_pfIn); if (m_pfOut) fclose(m_pfOut);

    for (short i=0;i<24;i++) {
        if (m_apcHsRef[i]) {
           delete [] m_apcHsRef[i];
        }
        
        if (m_apfHsRef[i]) fclose(m_apfHsRef[i]);
    }

    delete m_apcHsRef;
}


void GetHRefSeq::run_SV(const char *pcInput, const char *pcRefFilePath,
                     const char *pcOutput, int nSize)
{
    fprintf(stdout,"Open files ...\n");
    if (!OpenFiles(pcInput,pcRefFilePath, pcOutput)) return;

    fprintf(stdout,"Load Human Reference Sequence ...\n");
    LoadHRefSeq();

    fprintf(stdout,"Process Input File...\n");

    char acbuf[5096],acOut[5096]; char *pChr=NULL;
    short sChroIdx1, sChroIdx2; int nOffset1, nOffset2, nRange1, nRange2, nVal;
    char *lpcLeftSeq[3],*lpcRightSeq[3];

#if _TESTING
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
#endif

    while (!feof(m_pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),m_pfIn)) break;
        
	pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

        if (acbuf[0]==0) continue;
        if (acbuf[0]=='#'){
            fprintf(m_pfOut,"%s\tLeft_1\tLeft_2\tLeft_3\tRight_1\tRight_2\tRight_3\n",acbuf);
            continue;
        }
/*
#if _TESTING
        ulTest1++;
        if (ulTest1 == 1)//0000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = 36.3=== %u\n", ulTest2);
            ulTest1=0;
        }
#endif
*/
        strcpy(acOut,acbuf); pChr = strtok(acbuf,"\t");

        if (pChr[0] == 'X')
           sChroIdx1 = 22;
        else if (pChr[0] == 'Y')
           sChroIdx1 = 23;
        else
           sChroIdx1 = atoi(pChr)-1;

        pChr = strtok(NULL,"\t"); nOffset1 = atoi(pChr)-1;
        pChr = strtok(NULL,"\t"); nRange1 = atoi(pChr);
        pChr = strtok(NULL,"\t"); 

        pChr = strtok(NULL,"\t");

        if (pChr[0] == 'X')
           sChroIdx2 = 22;
        else if (pChr[0] == 'Y')
           sChroIdx2 = 23;
        else
           sChroIdx2 = atoi(pChr)-1;

        pChr = strtok(NULL,"\t"); nOffset2 = atoi(pChr)-1;
        pChr = strtok(NULL,"\t"); nRange2 = atoi(pChr);

        lpcLeftSeq[0] = new char[nRange1+1]; 
        memcpy(lpcLeftSeq[0],&m_apcHsRef[sChroIdx1][nOffset1],nRange1);
        lpcLeftSeq[0][nRange1]=0; 
        
        lpcLeftSeq[1] = new char[nSize+1]; 
        memcpy(lpcLeftSeq[1],&m_apcHsRef[sChroIdx1][nOffset1],nSize);        
        lpcLeftSeq[1][nSize]=0;

        lpcLeftSeq[2] = new char[nSize+1]; memset(lpcLeftSeq[2],0,nSize+1);
        nOffset1= (nOffset1+nRange1)-nSize;

        if (nOffset1 < 0){
            nVal = /*nSize - nRange1*/ nSize+nOffset1; nOffset1=0;}
        else
            nVal = nSize;

        memcpy(lpcLeftSeq[2],&m_apcHsRef[sChroIdx1][nOffset1],nVal);
        lpcLeftSeq[2][nVal]=0;  
        
        lpcRightSeq[0] = new char[nRange2+1];
        memcpy(lpcRightSeq[0],&m_apcHsRef[sChroIdx2][nOffset2],nRange2);        
        lpcRightSeq[0][nRange2]=0;

        lpcRightSeq[1] = new char[nSize+1];
        memcpy(lpcRightSeq[1],&m_apcHsRef[sChroIdx2][nOffset2],nSize);
        lpcRightSeq[1][nSize]=0;
        
        lpcRightSeq[2] = new char[nSize+1];
        nOffset2= (nOffset2+nRange2)-nSize;
        if (nOffset2 < 0){
            nVal = /*nSize - nRange2*/ nSize+nOffset2; nOffset2=0;}
        else
            nVal = nSize;

        memcpy(lpcRightSeq[2],&m_apcHsRef[sChroIdx2][nOffset2],nVal);
        lpcRightSeq[2][nVal]=0; 

       fprintf(m_pfOut,"%s\t%s\t%s\t%s\t%s\t%s\t%s\n",acOut,lpcLeftSeq[0],
                lpcLeftSeq[1],lpcLeftSeq[2],lpcRightSeq[0],lpcRightSeq[1],
                lpcRightSeq[2]);

       for (short i=0; i<3; i++){
           if (lpcLeftSeq[i]) delete[] lpcLeftSeq[i];
           if (lpcRightSeq[i])
               delete[] lpcRightSeq[i];
       }
    }

}


void GetHRefSeq::run_SnpIndel(const char *pcInput, const char *pcRefFilePath,
                              const char *pcOutput, int nSize)
{
    nSize++; //nSize + the current base refer by offset
    fprintf(stdout,"Open files ...\n");
    if (!OpenFiles(pcInput,pcRefFilePath, pcOutput)) return;

    fprintf(stdout,"Load Human Reference Sequence ...\n");
    LoadHRefSeq();
   
    fprintf(stdout,"Process Input File...\n");

    char acbuf[40960],acOut[40960]; char *pChr=NULL;
    short sChroIdx; int nOffset,nBuffSize = nSize;
    char *pcLeftSeq,*pcRightSeq;

//#if _TESTING
//    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
//#endif

    while (!feof(m_pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),m_pfIn)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

        if (acbuf[0]==0) continue;
        if (acbuf[0]=='C'){
            fprintf(m_pfOut,"%s\tLeft\tRight\n",acbuf);
            continue;
        }

/*#if _TESTING
        ulTest1++;
        if (ulTest1 == 1)//0000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u\n", ulTest2);
            ulTest1=0;
        }
#endif
*/

        strcpy(acOut,acbuf); pChr = strtok(acbuf,"\t");

        if (pChr[0] == 'X')
           sChroIdx = 22;
        else if (pChr[0] == 'Y')
           sChroIdx = 23;
        else
           sChroIdx = atoi(pChr)-1;

        pChr = strtok(NULL,"\t");
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr)-1;       

        pcRightSeq = new char[nBuffSize+1]; 

//if (ulTest2 == 1823)
//   fprintf(stdout,"nSize = %d\tm_anHRefSize=%d\tChroIdx=%d\tOffset = %d\n",nSize,m_anHsRefSize[sChroIdx],sChroIdx,nOffset);

        nSize = (m_anHsRefSize[sChroIdx]- nOffset)+1; 
        if (nSize > nBuffSize) nSize = nBuffSize;                    

        memcpy(pcRightSeq,&m_apcHsRef[sChroIdx][nOffset],nSize);
        pcRightSeq[nSize]=0;  

        pcLeftSeq = new char[nBuffSize+1]; 
        if (nSize != nBuffSize) nSize = nBuffSize; 

        nOffset = (nOffset-nSize)+1; if (nOffset < 0){nSize += nOffset; nOffset=0;}        
        memcpy(pcLeftSeq,&m_apcHsRef[sChroIdx][nOffset],nSize);
        pcLeftSeq[nSize]=0;                          
        
        fprintf(m_pfOut,"%s\t%s\t%s\n",acOut,pcLeftSeq,pcRightSeq);

        if (pcLeftSeq) delete[] pcLeftSeq;
        if (pcRightSeq) delete[] pcRightSeq;
        
    }
}


bool GetHRefSeq::OpenFiles(const char *pcInput, const char *pcRefFilePath,
                           const char *pcOutput)
{
    m_pfIn = fopen(pcInput,"r");

    if (!m_pfIn) {
        fprintf(stdout,"Failed to open %s ....",pcInput); return false;
    }
    
    m_pfOut = fopen(pcOutput,"w");
            
    if (!m_pfOut) {
        fclose(m_pfIn); 
        fprintf(stdout,"Failed to open %s ....",pcOutput); return false; 
    }

    char acFileName[1024];

    for (short i=0;i<24;i++)
    {
        if (i < 22){
            sprintf(acFileName,"%s/hs_ref_chr%d.fa",pcRefFilePath,i+1);
         }
        else if(i == 22)
            sprintf(acFileName,"%s/hs_ref_chrX.fa",pcRefFilePath);
        else
            sprintf(acFileName,"%s/hs_ref_chrY.fa",pcRefFilePath);

        m_apfHsRef[i] = fopen(acFileName,"r");

        if (m_apfHsRef[i]){
           fprintf(stdout,"Touch %s ...\n",acFileName);
        }
        else{
            fprintf(stdout,"Failed to open %s ....\n",acFileName);
            CloseAllFiles(); return false;
        }
    }
    return true;
}


inline void GetHRefSeq::CloseRefFiles()
{
    for (short i=0;i<24;i++) {
        if (m_apfHsRef[i]){ fclose(m_apfHsRef[i]); m_apfHsRef[i]=NULL;}
    }
}

void GetHRefSeq::CloseAllFiles()
{
    if (m_pfIn) {fclose(m_pfIn); m_pfIn=NULL;}
    if (m_pfOut) { fclose(m_pfOut); m_pfOut=NULL;}

    CloseRefFiles();  
}


void GetHRefSeq::LoadHRefSeq()
{
    struct stat st; size_t file_size;
    char acbuf[130]; char *pChr=NULL;
    int nLen,nTotalRecs; unsigned int unIdx;

//#if _TESTING
//    for (short i=0;i<1;i++)
//#else
    for (short i=0;i<24;i++)
//#endif
    {
        fstat(fileno(m_apfHsRef[i]), &st);
        file_size = (size_t)st.st_size;
        fgets(acbuf, sizeof(acbuf),m_apfHsRef[i]);

        nLen = strlen(acbuf);
        nTotalRecs = file_size - nLen;

        m_apcHsRef[i] = new char[nTotalRecs];
        unIdx=0;

        while(!feof(m_apfHsRef[i]))
        {
            if (!fgets(acbuf, sizeof(acbuf),m_apfHsRef[i])) break;

            pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
            pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

            nLen = strlen(acbuf); memcpy(m_apcHsRef[i]+unIdx,acbuf,nLen); unIdx+=nLen;
        }
        m_apcHsRef[i][unIdx]=0; m_anHsRefSize[i]=unIdx;

//#if _TESTING
       fprintf(stdout,"Human Ref.File %d => ",i+1);
       for (int j=0; j < 70; j++)
       {
            fprintf(stdout,"%c",m_apcHsRef[i][j]);
       }
       fprintf(stdout,"...\n");
//#endif
   }

   CloseRefFiles();
}

