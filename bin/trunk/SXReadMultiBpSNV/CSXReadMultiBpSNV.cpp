#include "CSXReadMultiBpSNV.h"
#include "SXDBSNPgChecker.h"
#include "SXDBIndelgChecker.h"

const char *VERSION="1.0.0.0";

CSXReadMultiBpSNV::CSXReadMultiBpSNV()
{
    m_pcdbSnp=m_pcdbIndel=m_pcDensFilePath=NULL;
    m_nQScoreOfs=33; //Default to 33
}


CSXReadMultiBpSNV::~CSXReadMultiBpSNV()
{
    ClrMemoryMap();
}


void CSXReadMultiBpSNV::ClrMemoryMap()
{
     if (m_stSNVFile.fd > 0 && m_stSNVFile.pstSNV != (void*)-1){
        munmap(m_stSNVFile.pstSNV,m_stSNVFile.ulFileSize); close(m_stSNVFile.fd);
     }

     m_stSNVFile.fd = 0;
}


bool CSXReadMultiBpSNV::MapSNVFile()
{    
    m_stSNVFile.fd = open(m_pcInFile,O_RDONLY,S_IRUSR);
    if (m_stSNVFile.fd == -1) {fprintf(stdout,"Failed to open %s ...",m_pcInFile); return false;} 

    struct stat st; fstat(m_stSNVFile.fd,&st);

    m_stSNVFile.pstSNV = (stSNV*)mmap(NULL,st.st_size,PROT_READ,MAP_PRIVATE,m_stSNVFile.fd,0);

    if (m_stSNVFile.pstSNV == (void*)-1){
        fprintf(stdout, "Error: mmap failed: %s\n", strerror(errno)); return false;
    }

    m_stSNVFile.ulFileSize = st.st_size;

    return true;  
}


void CSXReadMultiBpSNV::PrintSNVFileContents(char *pcOutput)
{
    FILE *pfOut = fopen(pcOutput,"w");
    if (!pfOut){fprintf(stdout,"Failed to open %s ...",pcOutput); return;} 

    if (!MapSNVFile()) return; 

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    for (unsigned long i=0; i < (m_stSNVFile.ulFileSize/sizeof(stSNV)); i++)
    {         
         ulTest1++;

         if (ulTest1 == 1000000)
         {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
         }      

         fprintf(pfOut,"%u\t%u\t%c\t%c\t%u\t%u\t%u\t%c\t%c\n",
                 abs(m_stSNVFile.pstSNV[i].cChromosome),ComputeOffset(m_stSNVFile.pstSNV[i].acOffset),
                 m_stSNVFile.pstSNV[i].cRefBase,m_stSNVFile.pstSNV[i].cVarBase,m_stSNVFile.pstSNV[i].ucQryPos,
                 m_stSNVFile.pstSNV[i].ucQltyScore,m_stSNVFile.pstSNV[i].usPEnd,
                 (m_stSNVFile.pstSNV[i].usPEnd & 8192)==8192?'R':'L',
                 m_stSNVFile.pstSNV[i].cStrand == '0'?'-':'+');
    }

    fclose(pfOut);
}


void CSXReadMultiBpSNV::setMbpSNVListInput(const char *pcdbSnp, const char *pcdbIndel,
                                   const char *pcDensFilePath,
                                   unsigned int unMinSupportingReads,
                                   unsigned int unMaxSupportingReads,
                                   float fMinReadStrength,int nMinQScore,
                                   const char *pcSampleID,const char *pcGenomeVer)
{
    m_pcdbSnp=pcdbSnp; m_pcdbIndel=pcdbIndel; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength = fMinReadStrength; m_nMinQScore=nMinQScore;
    m_pcSampleID=pcSampleID; m_pcGenomeVer=pcGenomeVer;
    m_bMultiBpSNVList=true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern "C" int CompareSNVRecs(const void *a, const void *b)
{
    //if (((const stSNVRec *) a)->cChromosome != ((const stSNVRec *) b)->cChromosome)
    //   return ((const stSNVRec *) a)->cChromosome - ((const stSNVRec *) b)->cChromosome;

    if (((const stSNVRec *) a)->unOffset != ((const stSNVRec *) b)->unOffset)
        return ((const stSNVRec *) a)->unOffset - ((const stSNVRec *) b)->unOffset;

    if (((const stSNVRec *) a)->cRefBase != ((const stSNVRec *) b)->cRefBase)
        return ((const stSNVRec *) a)->cRefBase - ((const stSNVRec *) b)->cRefBase;

    return ((const stSNVRec *) a)->cVarBase - ((const stSNVRec *) b)->cVarBase;
}


void *SortSNVData(void *ptr)
{
    stSNVRecList *pSNVRecList = (stSNVRecList*)ptr;
    char tmpChr;

    for (int i=0; i < pSNVRecList->nUsed; i++)
    {
        //pSNVRecList->pSNVRec[i].unOffset = ComputeOffset(pSNVRecList->pSNVRec[i].acOffset);
        tmpChr = toupper(pSNVRecList->pSNVRec[i].cRefBase);
        pSNVRecList->pSNVRec[i].cRefBase = tmpChr;
        tmpChr = toupper(pSNVRecList->pSNVRec[i].cVarBase);
        pSNVRecList->pSNVRec[i].cVarBase = tmpChr;
    }
    qsort(pSNVRecList->pSNVRec,pSNVRecList->nUsed,sizeof(stSNVRec), CompareSNVRecs);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CSXReadMultiBpSNV::AllocMergeListSize()
{
    if (m_pSNVClustered == NULL) {
        merge_list_alloc = INIT_MERGED_LIST_SIZE;
        merge_list_used = 0;
        try {
                 m_pSNVClustered = (stSNVClusterRec*)malloc(sizeof(stSNVClusterRec)*merge_list_alloc);
            }
       catch (...) {
                 fprintf(stderr, "Error: insufficient memory for Pair Info list!\n");
                 fprintf(stderr, "Program terminated.\n");
                 exit(1);
       }

       return;
    }

    merge_list_alloc += EXT_MERGED_LIST_SIZE;
    try {
              m_pSNVClustered = (stSNVClusterRec*)realloc(m_pSNVClustered,sizeof(stSNVClusterRec)*merge_list_alloc);
    }
    catch (...) {
              fprintf(stderr, "Error: insufficient memory for Pair Info list\n");
              fprintf(stderr, "Program terminated.\n");
              exit(1);
        }
}


void CSXReadMultiBpSNV::SetPEndnStrandInfo(stSNVClusterRec &SNVClusterRec,unsigned short &usPEnd,char &cStrand)
{
    if (usPEnd != 0) SNVClusterRec.unPEndCnt++;

    unsigned char unChar='0';
    if ((usPEnd & 0x2000)==0x2000) unChar='1'; //Right hand side
     
    if (cStrand==0)  SNVClusterRec.unRvsCnt++;            
    else if (cStrand==1){ SNVClusterRec.unFwdCnt++; unChar |= 0x2;}

    SNVClusterRec.Add(&SNVClusterRec.pucMisc,unChar);
}


void CSXReadMultiBpSNV::SetMiscInfo(stSNVClusterRec &SNVClusterRec,
                                    unsigned short &usPEnd, char &cStrand, 
                                    unsigned short &usQryPos)
{
    if (usPEnd != 0) SNVClusterRec.unPEndCnt++;

    unsigned char ucSide=0, ucStrand=0;
    ucSide =  ((usPEnd & 0x2000)==0x2000)? 'R':'L'; //Right hand side
     
    if (cStrand==0)  {SNVClusterRec.unRvsCnt++; ucStrand = '-';}            
    else if (cStrand==1){ SNVClusterRec.unFwdCnt++; ucStrand = '+';}

    SNVClusterRec.AddMisc(ucSide,ucStrand,usQryPos);
}       


void CSXReadMultiBpSNV::GenerateSNVTbl()
{
    unsigned long ulTest1=0, ulTest2=0;

    fprintf(stdout,"Loading, converting and splitting...\n"); 
   
    int nChro; stSNVRecList aSNVRecList[24]; stSNVRecList *pSNVRecList=NULL;

    FILE *pfVRI = fopen("IndenticalVRBase.lst","w"); 

    m_pSNV = m_stSNVFile.pstSNV;   
                                 
    //for (unsigned long i=0; i < 1000000; i++) //Testing
    for (unsigned long i=0; i < (m_stSNVFile.ulFileSize/sizeof(stSNV)); i++) 
    {
        //printf("RefBase=%c\nVarBase=%c\n",m_pSNV[i].cRefBase,m_pSNV[i].cVarBase); return;
        if (m_pSNV[i].cRefBase == m_pSNV[i].cVarBase) {fprintf(pfVRI,"%d\n",i); continue;}

        nChro = abs(m_pSNV[i].cChromosome); pSNVRecList = &aSNVRecList[nChro-1];
       
        if (pSNVRecList->nUsed == pSNVRecList->nTotal){
            pSNVRecList->nTotal += SNVREC_SIZE;
            pSNVRecList->pSNVRec = (stSNVRec*)realloc(pSNVRecList->pSNVRec, sizeof(stSNVRec)*pSNVRecList->nTotal);
        }

        //memcpy(pSNVRecList->pSNVRec[pSNVRecList->nUsed].acOffset,m_pSNV[i].acOffset,4);
        memcpy(&pSNVRecList->pSNVRec[pSNVRecList->nUsed].unOffset,m_pSNV[i].acOffset,4);
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cRefBase = m_pSNV[i].cRefBase;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cVarBase = m_pSNV[i].cVarBase;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQryPos = m_pSNV[i].ucQryPos;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQltyScore = m_pSNV[i].ucQltyScore;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].usPEnd = m_pSNV[i].usPEnd;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cStrand = m_pSNV[i].cStrand;
        pSNVRecList->nUsed++;
    }  
   
    fprintf(stdout,"Total SNV Recs = %u\n",m_stSNVFile.ulFileSize/sizeof(stSNV));       
    ClrMemoryMap(); fclose(pfVRI);
     
    fprintf(stdout,"sorting...\n");
    //return;

    for (int i=0; i<24; i++){
         if (pthread_create(&(aSNVRecList[i].ThreadID),NULL,SortSNVData,(void*)(&aSNVRecList[i])) != 0){
             fprintf(stdout,"Failed to spawn SortSNVData Thread %d ...\n",i);
             exit(0);
         }
    }

    for (int i=0; i<24; i++){
         if (pthread_join(aSNVRecList[i].ThreadID,NULL)!=0){
             fprintf(stdout,"Failed to perform Thread Join %d ...\n",i);
             exit(0);
         }
    }

    fprintf(stdout,"alloc memory for clustering...\n");

    AllocMergeListSize();

    fprintf(stdout,"clustering...\n");
    
    char cChromosome_Curr,cRefBase_Curr,cVarBase_Curr;
    char cRef_ACGT,cVar_ACGT,cStrand;
    unsigned int unOffset_Curr,unOffset;
    unsigned short usQP, usQS; unsigned short usPEnd=0;

    pSNVRecList=&aSNVRecList[0];

    unsigned short ausPEType[4]={49152,32768,16384,0};

    for (int idx=0; idx < 4; idx++)
    {
       if ((pSNVRecList->pSNVRec[0].usPEnd & ausPEType[idx])==ausPEType[idx]){
          usPEnd = pSNVRecList->pSNVRec[0].usPEnd - ausPEType[idx];
          break;
       }
    }

    usQP = pSNVRecList->pSNVRec[0].ucQryPos;
    usQS = pSNVRecList->pSNVRec[0].ucQltyScore-m_nQScoreOfs;
    cChromosome_Curr = 1;
    cStrand = pSNVRecList->pSNVRec[0].cStrand;
    unOffset_Curr = pSNVRecList->pSNVRec[0].unOffset;
    cRefBase_Curr = pSNVRecList->pSNVRec[0].cRefBase;
    cVarBase_Curr = pSNVRecList->pSNVRec[0].cVarBase;

    m_pSNVClustered[merge_list_used].InitRec();
    m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
    m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
    m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
    m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
    m_pSNVClustered[merge_list_used].ulTotalQP = usQP;
    m_pSNVClustered[merge_list_used].ulTotalQS = usQS;
    m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
    m_pSNVClustered[merge_list_used].unCnt = 1;
    //m_pSNVClustered[merge_list_used].bKeep = true;

    SetMiscInfo(m_pSNVClustered[merge_list_used],usPEnd,cStrand,usQP);

    if ((cRefBase_Curr!='A' && cRefBase_Curr!='C' && cRefBase_Curr!='G' && cRefBase_Curr!='T' &&
         cRefBase_Curr!='-')|| (cVarBase_Curr!='A' && cVarBase_Curr!='C' && cVarBase_Curr!='G' &&
         cVarBase_Curr!='T' && cVarBase_Curr!='-'))
    {
        //m_pSNVClustered[merge_list_used].bKeep=false;
        //LogInvalidRec(pfLog,pTmp[0]);
    }

   for (int i=0; i < 24; i++)
   {
       pSNVRecList=&aSNVRecList[i];

       for (int j=0; j < pSNVRecList->nUsed; j++)
       {

           if (i==0 && j==0) continue;

           if ((pSNVRecList->pSNVRec[j].cRefBase!='A' && pSNVRecList->pSNVRec[j].cRefBase!='C'
                && pSNVRecList->pSNVRec[j].cRefBase!='G'&& pSNVRecList->pSNVRec[j].cRefBase!='T'
                && pSNVRecList->pSNVRec[j].cRefBase!='-')||
               (pSNVRecList->pSNVRec[j].cVarBase!='A' && pSNVRecList->pSNVRec[j].cVarBase!='C'
                && pSNVRecList->pSNVRec[j].cVarBase!='G' && pSNVRecList->pSNVRec[j].cVarBase!='T'
                && pSNVRecList->pSNVRec[j].cVarBase!='-')
              )
           {
              // LogInvalidRec(pfLog,pTmp[i]); 
              continue;
           }

           cRef_ACGT = pSNVRecList->pSNVRec[j].cRefBase;
           cVar_ACGT = pSNVRecList->pSNVRec[j].cVarBase;
           unOffset = pSNVRecList->pSNVRec[j].unOffset;
           usQP = pSNVRecList->pSNVRec[j].ucQryPos;
           usQS = pSNVRecList->pSNVRec[j].ucQltyScore-m_nQScoreOfs;

           for (int idx=0; idx < 4; idx++)
           {
               if ((pSNVRecList->pSNVRec[j].usPEnd & ausPEType[idx])==ausPEType[idx]) {
                  usPEnd = pSNVRecList->pSNVRec[j].usPEnd - ausPEType[idx]; break;
               }
           }

           cStrand = pSNVRecList->pSNVRec[j].cStrand;

           if (cChromosome_Curr != i+1||unOffset_Curr != unOffset||
               cRefBase_Curr != cRef_ACGT||cVarBase_Curr != cVar_ACGT)
           {
               cChromosome_Curr = i+1; unOffset_Curr = unOffset;
               cRefBase_Curr = cRef_ACGT; cVarBase_Curr = cVar_ACGT;

               merge_list_used++; 
               if (merge_list_used == merge_list_alloc) AllocMergeListSize();

               m_pSNVClustered[merge_list_used].InitRec();
               m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
               m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
               m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
               m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
               m_pSNVClustered[merge_list_used].ulTotalQP = usQP;
               m_pSNVClustered[merge_list_used].ulTotalQS = usQS;
               m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
               //m_pSNVClustered[merge_list_used].bKeep = true;
          }
          else
          {
               m_pSNVClustered[merge_list_used].ulTotalQP += usQP;
               m_pSNVClustered[merge_list_used].ulTotalQS += usQS;

               if (m_pSNVClustered[merge_list_used].ucMaxQS < usQS)
                   m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
          }

          m_pSNVClustered[merge_list_used].unCnt++;
          SetMiscInfo(m_pSNVClustered[merge_list_used],usPEnd,cStrand,usQP);

          ulTest1++;
          if (ulTest1 == 1000000)
          {
              ulTest2++;
              fprintf(stdout,"Total Recs gone through clustering = %uM\n", ulTest2);
              ulTest1=0;
          }
       }
    }

    fprintf(stdout,"Total Clustered Recs = %u\n", merge_list_used);
    return;
}


void CSXReadMultiBpSNV::Cluster()
{
    if (!MapSNVFile()) return;

    FILE *pfdbSnp=NULL, *pfdbIndel=NULL;

    char acFile[strlen(m_pcDensFilePath)+16];
    m_apfRD = new FILE*[24]; for (int i=0; i<24; i++){m_apfRD[i]=NULL;}
    pfdbSnp = fopen(m_pcdbSnp,"r");if (!pfdbSnp)
    {
        fprintf(stderr,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;
    }

    pfdbIndel = fopen(m_pcdbIndel,"r");if (!pfdbSnp)
    {
        fprintf(stderr,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;
    }

    for (int i=0; i<24;i++)
    {
        sprintf(acFile,"%s/Chr%02d_read_den",m_pcDensFilePath,i+1);
        m_apfRD[i] = fopen(acFile,"r");
        if (!m_apfRD[i]){
           fprintf(stderr,"Failed to open %s ...\n", acFile); goto Exit;
        }
        else
           fprintf(stderr,"Touch %s ...\n",acFile);
    }

    fprintf(stdout,"Generating dbSNP Map...\n");
    if (!GenerateDBSNPMap(pfdbSnp)) goto Exit;
    fprintf(stdout,"Generating dbINDEL Map...\n");
    if (!GenerateDBINDELMap(pfdbIndel)) goto Exit;

    fprintf(stdout,"Generating SNV Table...\n"); GenerateSNVTbl();
    fprintf(stdout,"Printing SNP, INS, DEL Lists...\n"); PrintMultiBpSNVList();

Exit:     
    
    if (pfdbSnp) fclose(pfdbSnp); if (pfdbIndel) fclose(pfdbIndel);
    ClrDBSNPMap(); ClrDBINDELMap();

    for (int j=0;j<24; j++){
        if (m_apfRD[j]) fclose(m_apfRD[j]);
    }

    delete[] m_apfRD; 
}


char g_acMonth[12][4]= {"Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec"};

inline void CSXReadMultiBpSNV::PrintRptDateTime(FILE *pf)
{   
    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);
    
    fprintf(pf,"%02d-%s-20%02d; %02d:%02d:%02d\n",
               ptm->tm_mday,g_acMonth[ptm->tm_mon],ptm->tm_year-100,
               ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
}


inline void CSXReadMultiBpSNV::PrintRemarksParam(FILE *pfParam)
{
    FILE *pf = fopen("Param.ini","r");
    if (!pf) return;

    char acbuf[1024]; size_t unSize=sizeof(acbuf);
    char *pChr=NULL;

    while(!feof(pf))
    {
        if (!fgets(acbuf,unSize,pf)){break;}
        if (acbuf[0]=='#') continue;

        pChr = strtok(acbuf,"=");

        if (strcmp(pChr,"FRemarks")==0)
        {
            pChr = strtok(NULL,"\n");

            if (pChr)
                fprintf(pfParam,"%s\n",pChr);
        }
    }

    fprintf(pfParam,"\n");
    if (pf) {fclose(pf);}
}


inline double CSXReadMultiBpSNV::CalcAvgQS(unsigned long ulTotalQS, unsigned long ulTotalCnt)
{
    if (ulTotalCnt == 0) return 0;
    return (ulTotalQS == 0)? 0:(double)ulTotalQS/ulTotalCnt;
}


inline void CSXReadMultiBpSNV::OutputChrom(FILE *pf, unsigned int unChrom)
{
    if (unChrom < 23)
        fprintf(pf,"%d",unChrom);
    else if (unChrom == 23)
        fprintf(pf,"X");
    else
        fprintf(pf,"Y");
}


unsigned char g_acBuf[2];unsigned int g_unReads;
inline unsigned int CSXReadMultiBpSNV::GetReadDensity(FILE *pf, unsigned int uioffset)
{
    if (!pf) return 0;
    uioffset--;
    fseek(pf,uioffset*2,SEEK_SET);
    fread(g_acBuf,2,1,pf);

    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|= g_acBuf[0];
    
    return g_unReads;
}
    

inline unsigned int CSXReadMultiBpSNV::GetReadDensityEx(FILE *pf, unsigned int uioffset)
{
    if (!pf) return 0;

    unsigned int unRead=0, unTotalReads=0;
    uioffset--;

    for (int i=1;i<2;i++)
    {
         fseek(pf,(uioffset-i)*2,SEEK_SET);
         fread(g_acBuf,2,1,pf);

         unRead = g_acBuf[1];unRead<<=8;
         unRead|= g_acBuf[0];
         unTotalReads += unRead;
    }

    for (int i=1;i<2;i++)
    {
         fseek(pf,(uioffset+i)*2,SEEK_SET);
         fread(g_acBuf,2,1,pf);

         unRead = g_acBuf[1];unRead<<=8;
         unRead|= g_acBuf[0];
         unTotalReads += unRead;
    }

    g_unReads = (unsigned int)((double)((double)(unTotalReads)/2)+0.5);

    return g_unReads;
}


inline double CSXReadMultiBpSNV::Percentage(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator*100;
}


inline double CSXReadMultiBpSNV::Division(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator;
}


void CSXReadMultiBpSNV::PrintMultiBpSNVList()
{
    /////////////////////////////

    stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astDELStats[m_unMaxSupportingReads+1];
    int nSNPIdx,nINSIdx,nDELIdx;

    unsigned long ulTest1=0, ulTest2=0;
    unsigned int unChroIdx,unReadDensity=0;

    bool bdbSNP=false; char acRsid[200];

    float fAvgQryPos,fAvgQScore;

    FILE *pfSNP,*pfINS,*pfDEL,*pfSNPStats,*pfINSStats,*pfDELStats;
    pfSNP=pfINS=pfDEL=pfSNPStats=pfINSStats=pfDELStats=0;

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[24+strlen(m_pcSampleID)];
    sprintf(acFileName,"%s_mbsnp_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}


    sprintf(acFileName,"%s_mbsnp_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    sprintf(acFileName,"%s_mbins_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_mbins_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////
    sprintf(acFileName,"%s_mbdel_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_mbdel_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    fprintf(pfSNP,"#REPORT NAME\tSNP Distribution\n");
    fprintf(pfSNP,"#PROJECT NAME\n");
    fprintf(pfSNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfSNP,"#LANE NO\tALL\n");
    fprintf(pfSNP,"#GENERATED AT\t"); PrintRptDateTime(pfSNP);
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfSNP,"#REMARKS\t");PrintRemarksParam(pfSNP);
    fprintf(pfSNP,"#FILTER\t>=1\n\n\n");
    fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"Chromosome\tOffset\tNucleotide_Variant\tFwd_SNP_Reads"
                  "\tRvs_SNP_Reads\tTotal_SNP_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\tMisc.\n");

    fprintf(pfINS,"#REPORT NAME\tINS Distribution\n");
    fprintf(pfINS,"#PROJECT NAME\n");
    fprintf(pfINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfINS,"#LANE NO\tALL\n");
    fprintf(pfINS,"#GENERATED AT\t"); PrintRptDateTime(pfINS);
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfINS,"#REMARKS\t");PrintRemarksParam(pfINS);
    fprintf(pfINS,"#FILTER\t>=1\n\n\n");

    fprintf(pfINS,"Chromosome\tOffset\tInserted_Base\tFwd_INS_Reads"
                  "\tRvs_INS_Reads\tTotal_INS_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\tMisc.\n");

    fprintf(pfDEL,"#REPORT NAME\tDEL Distribution\n");
    fprintf(pfDEL,"#PROJECT NAME\n");
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfDEL,"#REMARKS\t");PrintRemarksParam(pfDEL);
    fprintf(pfDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfDEL,"Chromosome\tOffset\tDeleted_Base\tFwd_DEL_Reads"
                  "\tRvs_DEL_Reads\tTotal_DEL_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\tMisc.\n");

    for (unsigned long i=0; i < merge_list_used; i++)
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Clustered Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }

        if (m_pSNVClustered[i].cChromosome == 0) continue;

        if (m_pSNVClustered[i].cRefBase!='-' && m_pSNVClustered[i].cVarBase!='-')
        {
            bdbSNP=OutputSNPNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                  m_pSNVClustered[i].cVarBase, pfSNP,acRsid);

            //if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            //if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            //if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfSNP,m_pSNVClustered[i].cChromosome);
            fprintf(pfSNP,"\t%u\t%c>%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,
                    m_pSNVClustered[i].unCnt);

            fprintf(pfSNP,"\t%u\t%u\t%.2f\t%.2f\t%s\t%s\n",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid,m_pSNVClustered[i].pucMisc);

            nSNPIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astSNPStats[nSNPIdx].unCnt++;
            if (bdbSNP) astSNPStats[nSNPIdx].undbSNP++;
            astSNPStats[nSNPIdx].ulTotalReadDens += unReadDensity;
  
        }
        else if (m_pSNVClustered[i].cRefBase =='-' )
        {
            bdbSNP=OutputINDELNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                    m_pSNVClustered[i].cVarBase, pfINS,acRsid);

            //if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset);

            //if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            //if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfINS,m_pSNVClustered[i].cChromosome);

            fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,
                    m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,m_pSNVClustered[i].unCnt);

            fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s\t%s\n",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid,m_pSNVClustered[i].pucMisc);

            nINSIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astINSStats[nINSIdx].unCnt++;
            if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
            astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;
        }
        else
        {
            bdbSNP=OutputINDELNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                    m_pSNVClustered[i].cRefBase,pfDEL,acRsid);

            //if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            //if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            //if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfDEL,m_pSNVClustered[i].cChromosome);
            fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,m_pSNVClustered[i].unCnt);

            fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s\t%s\n",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid,m_pSNVClustered[i].pucMisc);

            nDELIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astDELStats[nDELIdx].unCnt++;
            if (bdbSNP) astDELStats[nDELIdx].undbSNP++;
            astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;
        }

        free(m_pSNVClustered[i].pucMisc);  m_pSNVClustered[i].pucMisc = NULL;
    }


    for (unsigned long i=0; i < merge_list_used; i++) m_pSNVClustered[i].ClrRec();


    fprintf(stdout,"Printing AvgQScore Statistic Report...\n");

    if (astSNPStats[0].unCnt > 0) PrintAvgQSStatsRpt(pfSNPStats,astSNPStats,"SNP");
    if (astINSStats[0].unCnt > 0) PrintAvgQSStatsRpt(pfINSStats,astINSStats,"INS");
    if (astDELStats[0].unCnt > 0) PrintAvgQSStatsRpt(pfDELStats,astDELStats,"DEL");    

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfINS) fclose(pfINS);
    if (pfDEL) fclose(pfDEL); if (pfSNPStats) fclose(pfSNPStats);
    if (pfINSStats) fclose(pfINSStats); if (pfDELStats) fclose(pfDELStats);
}


void CSXReadMultiBpSNV::PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{
    fprintf(pf,"#REPORT NAME\t%s Avg. QScore Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#REMARKS\t");PrintRemarksParam(pf);fprintf(pf,"\n");
    fprintf(pf,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
               m_unMinSupportingReads,m_unMaxSupportingReads,m_fMinReadStrength);


    fprintf(pf,"Supporting_Reads\t%s_Counts\tAvg_Read_Density\tdbSNP_Count\t"
               "dbSNP_%%\t\t\t"
               "Accum_SupportingReads\tAccum_%s_Counts\tAccum_%s_Counts_%%\t"
               "Accum_dbSNP_Count\tAccum_dbSNP_%%\n",pcType,pcType,pcType);


   m_stSNVAnalysis.ulTotalSNP=0; m_stSNVAnalysis.ulTotalSNP_dbsnp=0;

    unsigned long ulMax = m_unMaxSupportingReads+1;

    for (unsigned int i=0; i < ulMax; i++)
    {
        m_stSNVAnalysis.ulTotalSNP+=pstTblStats[i].unCnt;
        m_stSNVAnalysis.ulTotalSNP_dbsnp+=pstTblStats[i].undbSNP;
    }

    fprintf(pf,"=%u\t%u\t%.2f\t%u\t%.2f\t\t\t"
               ">=%u\t%u\t%.2f\t%u\t%.2f\n",
               m_unMinSupportingReads,pstTblStats[0].unCnt,
               Division(pstTblStats[0].ulTotalReadDens,pstTblStats[0].unCnt),
               pstTblStats[0].undbSNP,
               Percentage(pstTblStats[0].undbSNP,pstTblStats[0].unCnt),
               m_unMinSupportingReads,m_stSNVAnalysis.ulTotalSNP,
               (m_stSNVAnalysis.ulTotalSNP >0)?100.0:0.0,
               m_stSNVAnalysis.ulTotalSNP_dbsnp,
               Percentage(m_stSNVAnalysis.ulTotalSNP_dbsnp,m_stSNVAnalysis.ulTotalSNP));

    unsigned long ulAccumSnp,ulAccumdbSnp;

    ulAccumSnp=m_stSNVAnalysis.ulTotalSNP;
    ulAccumdbSnp=m_stSNVAnalysis.ulTotalSNP_dbsnp;

    for (unsigned int i=m_unMinSupportingReads; i < m_unMaxSupportingReads; i++){
        ulAccumSnp-=pstTblStats[i-1].unCnt;
        ulAccumdbSnp-=pstTblStats[i-1].undbSNP;

        fprintf(pf,"=%u\t%u\t%.2f\t%u\t%.2f\t\t\t"
                   ">=%u\t%u\t%.2f\t%u\t%.2f\n",
                   i+1,pstTblStats[i].unCnt,
                   Division(pstTblStats[i].ulTotalReadDens,pstTblStats[i].unCnt),
                   pstTblStats[i].undbSNP,
                   Percentage(pstTblStats[i].undbSNP,pstTblStats[i].unCnt),
                   i+1,ulAccumSnp,Percentage(ulAccumSnp,m_stSNVAnalysis.ulTotalSNP),
                   ulAccumdbSnp,Percentage(ulAccumdbSnp,ulAccumSnp));
    }


    fprintf(pf,">%u\t%u\t%.2f\t%u\t%.2f\t\t\t"
               ">=%u\t%u\t%.2f\t%u\t%.2f\n",
               m_unMaxSupportingReads,pstTblStats[m_unMaxSupportingReads].unCnt,
               Division(pstTblStats[m_unMaxSupportingReads].ulTotalReadDens,
                        pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].undbSNP,
               Percentage(pstTblStats[m_unMaxSupportingReads].undbSNP,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               m_unMaxSupportingReads+1,pstTblStats[m_unMaxSupportingReads].unCnt,
               Percentage(pstTblStats[m_unMaxSupportingReads].unCnt,
                          m_stSNVAnalysis.ulTotalSNP),
               pstTblStats[m_unMaxSupportingReads].undbSNP,
               Percentage(pstTblStats[m_unMaxSupportingReads].undbSNP,
                          pstTblStats[m_unMaxSupportingReads].unCnt));

    fprintf(pf,"Total\t%u\t\t%u\n",
            m_stSNVAnalysis.ulTotalSNP,m_stSNVAnalysis.ulTotalSNP_dbsnp);
    fprintf(pf,"\n\n\n");
}
