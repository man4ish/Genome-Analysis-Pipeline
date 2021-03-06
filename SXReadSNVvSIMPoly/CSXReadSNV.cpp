/* 
 * File:   CSXReadSNV.cc
 * Author: Manish
 * 
 * Created on November 17, 2009, 12:40 PM
 */

#include "CSXReadSNV.h"
#include "SXDBSNPgChecker.h"
#include "SXDBIndelgChecker.h"

#define _TESTING 0
#define VERSION "8.0.1.2"

/*
>gi|89161185|ref|NC_000001.9|NC_000001 Homo sapiens chromosome 1, reference assembly, complete sequence
>gi|89161199|ref|NC_000002.10|NC_000002 Homo sapiens chromosome 2, reference assembly, complete sequence
>gi|89161205|ref|NC_000003.10|NC_000003 Homo sapiens chromosome 3, reference assembly, complete sequence
>gi|89161207|ref|NC_000004.10|NC_000004 Homo sapiens chromosome 4, reference assembly, complete sequence
>gi|51511721|ref|NC_000005.8|NC_000005 Homo sapiens chromosome 5, reference assembly, complete sequence
>gi|89161210|ref|NC_000006.10|NC_000006 Homo sapiens chromosome 6, reference assembly, complete sequence
>gi|89161213|ref|NC_000007.12|NC_000007 Homo sapiens chromosome 7, reference assembly, complete sequence
>gi|51511724|ref|NC_000008.9|NC_000008 Homo sapiens chromosome 8, reference assembly, complete sequence
>gi|89161216|ref|NC_000009.10|NC_000009 Homo sapiens chromosome 9, reference assembly, complete sequence
>gi|89161187|ref|NC_000010.9|NC_000010 Homo sapiens chromosome 10, reference assembly, complete sequence
>gi|51511727|ref|NC_000011.8|NC_000011 Homo sapiens chromosome 11, reference assembly, complete sequence
>gi|89161190|ref|NC_000012.10|NC_000012 Homo sapiens chromosome 12, reference assembly, complete sequence
>gi|51511729|ref|NC_000013.9|NC_000013 Homo sapiens chromosome 13, reference assembly, complete sequence
>gi|51511730|ref|NC_000014.7|NC_000014 Homo sapiens chromosome 14, reference assembly, complete sequence
>gi|51511731|ref|NC_000015.8|NC_000015 Homo sapiens chromosome 15, reference assembly, complete sequence
>gi|51511732|ref|NC_000016.8|NC_000016 Homo sapiens chromosome 16, reference assembly, complete sequence
>gi|51511734|ref|NC_000017.9|NC_000017 Homo sapiens chromosome 17, reference assembly, complete sequence
>gi|51511735|ref|NC_000018.8|NC_000018 Homo sapiens chromosome 18, reference assembly, complete sequence
>gi|42406306|ref|NC_000019.8|NC_000019 Homo sapiens chromosome 19, reference assembly, complete sequence
>gi|51511747|ref|NC_000020.9|NC_000020 Homo sapiens chromosome 20, reference assembly, complete sequence
>gi|51511750|ref|NC_000021.7|NC_000021 Homo sapiens chromosome 21, reference assembly, complete sequence
>gi|89161203|ref|NC_000022.9|NC_000022 Homo sapiens chromosome 22, reference assembly, complete sequence
>gi|89161218|ref|NC_000023.9|NC_000023 Homo sapiens chromosome X, reference assembly, complete sequence
>gi|89161220|ref|NC_000024.8|NC_000024 Homo sapiens chromosome Y, reference assembly, complete sequence
 */

 unsigned int g_unGINums_G36[24] = {89161185,89161199,89161205,89161207,51511721,89161210,
                                    89161213,51511724,89161216,89161187,51511727,89161190,
                                    51511729,51511730,51511731,51511732,51511734,51511735,
                                    42406306,51511747,51511750,89161203,89161218,89161220};

 
 unsigned int g_unGINums_G37[24] = {224589800,224589811,224589815,224589816,224589817,224589818,
                                    224589819,224589820,224589821,224589801,224589802,224589803,
                                    224589804,224589805,224589806,224589807,224589808,224589809,
                                    224589810,224589812,224589813,224589814,224589822,224589823};

 extern char * g_geneIds[20];

 CSXReadSNV::CSXReadSNV(){
   m_pcInFile=m_pcOutFile=m_pcStatsFile=m_pcFreqFile=m_pcChroFreqFile=
   m_pcStatsFile=m_pcSNPStatsFile=m_pcINSStatsFile=m_pcDELStatsFile=
   m_pcChroStrandFile=m_pcChromOffsetFile=m_pcSNVTbl=m_pcNewFile=
   m_pcBinnedFile=m_pcAnnotateFile=m_pcdbSnp=m_pcdbIndel=m_pcSampleID=0;

   m_punGINums=0;

   m_bOut2Screen=m_bFQryPos=m_bFQltyScore=m_bMemMapped=m_bStateRpt=
   m_bSNVLst=m_bFilter=/*m_bSNPDensity=m_bINDELDensity=*/
   m_bDensity=m_bSNVTrace=m_bAvgQScoreLst=false; 

   m_nFQryPos1=m_nFQryPos2=m_nFQltyScore1=m_nFQltyScore2=-999;
   m_unMinSupportingReads=0; m_unMaxSupportingReads=70; m_unMaxReadDensity = 100;
   //m_fAvgQScore=
   m_fMinReadStrength=0;
   m_nQScoreOfs=33; //Default to 33   
   m_eSNV = eALL; m_pSXAnnotate = NULL;
   GetInitParam();

   #if _TESTING
        INIT_MERGED_LIST_SIZE = 500;
        EXT_MERGED_LIST_SIZE = 100;
   #else
        INIT_MERGED_LIST_SIZE = 500000000;
        EXT_MERGED_LIST_SIZE = 100000000;
   #endif        
}


CSXReadSNV::~CSXReadSNV() {
    if (m_pSXAnnotate) delete m_pSXAnnotate;
}


void CSXReadSNV::GetInitParam()
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

        if (strcmp(pChr,"QScoreOffset")==0)
        {
            pChr = strtok(NULL,"\n");

            if (pChr)
                m_nQScoreOfs = atoi(pChr);
        }
    }

    if (pf) {fclose(pf);}
}


void CSXReadSNV::setFilterParams(int nQryPos1,int nQryPos2,int nQltyScore1,
                                 int nQltyScore2,char cSNVType,const char* pcFile)
{
    //m_bFilter = true;
    m_nFQryPos1 = nQryPos1-1;
    m_nFQryPos2 = nQryPos2+1;
    m_nFQltyScore1 = nQltyScore1-1;
    m_nFQltyScore2 = nQltyScore2+1;
    m_pcNewFile = pcFile;

    m_eSNV = (eSNVType)cSNVType;
}


void CSXReadSNV::setSNVTraceInputs(unsigned int unSupportingReads,char cSNV)
{
    m_unMinSupportingReads = unSupportingReads;
    m_eSNV = (eSNVType)cSNV;
}


void CSXReadSNV::setQryPos(int nQryPos1, int nQryPos2, const char* pcFile)
{
    m_bFQryPos = true;
    m_nFQryPos1 = nQryPos1-1;
    m_nFQryPos2 = nQryPos2+1;
    m_pcNewFile = pcFile;
}



void CSXReadSNV::setQltyScore(int nQltyScore1,int nQltyScore2,const char* pcFile)
{
    m_bFQltyScore = true;
    m_nFQltyScore1 = nQltyScore1-1;
    m_nFQltyScore2 = nQltyScore2+1;
    m_pcNewFile = pcFile;
}


void CSXReadSNV::setAvgQSLstInputs(const char *pcdbSnp, const char *pcdbIndel,
                                   const char *pcDensFilePath, unsigned int unMinSupportingReads,
                                   unsigned int unMaxSupportingReads, float fMinReadStrength,
                                   int nMinQScore, const char *pcSampleID,const char *pcGenomeVer,
                                   int nQScoreOfs,const char *pcHsRefFilePath)
{
    m_pcdbSnp=pcdbSnp; m_pcdbIndel=pcdbIndel; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads; m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength=fMinReadStrength; m_nMinQScore=nMinQScore; m_pcSampleID=pcSampleID; 
    m_pcGenomeVer=pcGenomeVer; m_nQScoreOfs=nQScoreOfs; m_pcHsRefFilePath=pcHsRefFilePath;
    m_bAvgQScoreLst=true;
}


void CSXReadSNV::setSNVLstInputs(const char *pcBinnedFile, const char *pcExceptionFile,
                                 const char *pcAnnotateFile,
                                 unsigned int unMinSupportingReads,
                                 unsigned int unMaxSupportingReads,
                                 char cSNVType, const char *pcSampleID,
                                 const char *pcGenomeVer)
{
    m_pcBinnedFile=pcBinnedFile; m_pcExceptionFile=pcExceptionFile;
    m_pcAnnotateFile=pcAnnotateFile;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;    
    m_eSNV = (eSNVType)cSNVType;
    m_pcSampleID=pcSampleID;
    m_bSNVLst = true;
    m_pcGenomeVer=pcGenomeVer;
}


void CSXReadSNV::setDensityInput(const char* pcFile, char cSNVType, const char* pcSampleID)
{
    m_pcInFile=pcFile; m_eSNV=(eSNVType)cSNVType; 
    m_pcSampleID=pcSampleID; m_bDensity=true; 
}


void CSXReadSNV::run()
{   
    int fd =0;     
    if (m_bSNVLst) {OutputSNVList(); goto Exit;}

    fd=open(m_pcInFile, O_RDONLY, S_IRUSR);
    if (fd == -1) {fprintf(stderr, "File %s not found...\n",m_pcInFile); goto Exit;}

    struct stat st; fstat(fd, &st);

    m_file_size = (size_t)st.st_size; m_ulTotalRecs = m_file_size/sizeof(stSNV);

    m_pSNV=(stSNV*)mmap(NULL,m_file_size,PROT_READ,MAP_PRIVATE,fd,0);

    if (m_pSNV == (void *)-1) {
        fprintf(stderr, "Error: mmap failed: %s\n", strerror(errno));goto Exit;
    }

    m_bMemMapped = true;
   
    if (m_bFilter) {
        fprintf(stdout,"Filtering SNV List ....\n");

        if (m_eSNV == eSNP) FilterSNVList_SNP();
        else if (m_eSNV==eINS) FilterSNVList_INS();
        else if (m_eSNV==eDEL) FilterSNVList_DEL();
        else FilterSNVList();

        goto ClrMMap;
    }
    else if (m_bFilterx){
        fprintf(stdout,"Filtering SNV List ....\n");

        if (m_eSNV == eSNP) FilterSNVListx_SNP();
        else if (m_eSNV==eINS) FilterSNVListx_INS();
        else if (m_eSNV==eDEL) FilterSNVListx_DEL();
        else FilterSNVListx();

        goto ClrMMap;
    }

    if (m_bFQryPos) {fprintf(stderr,"Filtering SNV List By Query Position....\n");FilterSNVList2FileByQP(); goto ClrMMap;}
    if (m_bFQltyScore) {fprintf(stderr,"Filtering SNV List By Quality Score....\n");FilterSNVList2FileByQS(); goto ClrMMap;}

    if (m_pcStatsFile){fprintf(stderr,"Generating SNV Stats Report ....\n");OutputSNVStats('a',m_pcStatsFile);}
    if (m_pcSNPStatsFile){fprintf(stderr,"Generating SNP Stats Report ....\n");OutputSNVStats('s',m_pcSNPStatsFile);}
    if (m_pcINSStatsFile){fprintf(stderr,"Generating INS Stats Report ....\n");OutputSNVStats('i',m_pcINSStatsFile);}
    if (m_pcDELStatsFile){fprintf(stderr,"Generating DEL Stats Report ....\n");OutputSNVStats('d',m_pcDELStatsFile);}

    if (m_bAvgQScoreLst){OutputAvgQScoreTable(); goto Exit;}

    //if (m_bSNPDensity){OutputSNPDensity1K();}
    //if (m_bINDELDensity){OutputINDELDensity1K();}

    if (m_bDensity){OutputDensity();} 

    if (m_pcOutFile) {fprintf(stderr,"Generating SNV ascii file ....\n");OutputInVal2File();}
    if (m_bOut2Screen) OutputInVal2Scrn();

    if (m_pcFreqFile) {fprintf(stderr,"Generating SNV Freq Report ....\n");OutputSNVFreq();}
    if (m_pcChroFreqFile) {fprintf(stderr,"Generating Chromosome Freq Report ....\n");OutputChroFreq();}
    if (m_pcChroStrandFile) {fprintf(stderr,"Generating Chromosome Strand Freq Report ....\n");OutputChroStrand();}
    if (m_pcChromOffsetFile) {fprintf(stderr,"Generating Chromosome Offset Freq Report ....\n");OutputChromOffset();}

    if (m_bSNVTrace){fprintf(stderr,"Tracing Supporting data ....\n");TraceSNVSupportingReads();}

ClrMMap:
    ClrMemoryMap();

Exit:
    if (m_pSNVClustered) free(m_pSNVClustered);//delete[] m_pSNVClustered;
    if (fd) close(fd);
}


void CSXReadSNV::FilterSNVList()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
       bQScore=bQPos=false;

       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
           bQPos = true;

       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
           bQScore = true;

       if (bQScore&&bQPos)
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVList_SNP()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_pSNV[i].cRefBase == '-' || m_pSNV[i].cVarBase == '-') continue;
        if (m_pSNV[i].cChromosome == '\0') continue;

       bQScore=bQPos=false;

       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
           bQPos = true;

       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
           bQScore = true;

       if (bQScore&&bQPos)
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVList_INS()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
       if (m_pSNV[i].cRefBase != '-' ) continue;
       if (m_pSNV[i].cChromosome == '\0') continue;

       bQScore=bQPos=false;

       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
           bQPos = true;

       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
           bQScore = true;

       if (bQScore&&bQPos)
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}

void CSXReadSNV::FilterSNVList_DEL()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
       if (m_pSNV[i].cVarBase != '-' ) continue;
       if (m_pSNV[i].cChromosome == '\0') continue;

       bQScore=bQPos=false;

       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
           bQPos = true;

       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
           bQScore = true;

       if (bQScore&&bQPos)
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVListx()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_pSNV[i].cChromosome == '\0')
            continue;

       bQScore=bQPos=true;

       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
           bQPos = false;

       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
           bQScore = false;

       if (bQScore||bQPos)
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVListx_SNP()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_pSNV[i].cChromosome == '\0') continue;

        bQScore=bQPos=true;

        if (m_pSNV[i].cRefBase != '-' && m_pSNV[i].cVarBase != '-')
        {
            if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
                bQPos = false;

            nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

            if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
                bQScore = false;
        }
        if (bQScore||bQPos)
            fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVListx_INS()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_pSNV[i].cChromosome == '\0') continue;

        bQScore=bQPos=true;

        if (m_pSNV[i].cRefBase == '-')
        {
            if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
                bQPos = false;

            nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

            if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
                bQScore = false;
        }

        if (bQScore||bQPos)
            fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVListx_DEL()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    bool bQScore,bQPos; int nScore;

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_pSNV[i].cChromosome == '\0') continue;

        bQScore=bQPos=true;

        if (m_pSNV[i].cVarBase == '-'){
            if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
                bQPos = false;

            nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;

            if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
                bQScore = false;
        }
        if (bQScore||bQPos)
            fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
    }

    fclose(pf);
}


void CSXReadSNV::FilterSNVList2FileByQP()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {             
       if (m_pSNV[i].ucQryPos > m_nFQryPos1 && m_pSNV[i].ucQryPos < m_nFQryPos2)
       {
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
       }
    }
    fclose(pf);
}


void CSXReadSNV::FilterSNVList2FileByQS()
{
    FILE *pf=fopen(m_pcNewFile,"wb");
    if (!pf) {fprintf(stderr,"Can't open file %s...",m_pcNewFile); return;}

    int nScore;
    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
       nScore = m_pSNV[i].ucQltyScore-m_nQScoreOfs;
       if (nScore > m_nFQltyScore1 && nScore < m_nFQltyScore2)
       {
           fwrite(&m_pSNV[i],sizeof(stSNV),1,pf);
       }
    }
    fclose(pf);
}


void CSXReadSNV::ClrMemoryMap()
{
   if (m_bMemMapped)
   {
       int rc = munmap((char *)m_pSNV, m_file_size);
       if (rc == -1) {fprintf(stderr, "munmap failed...\n"); return;}
       m_bMemMapped = false;
   }
   else
   {
       if (m_pSNV) delete[] m_pSNV;
   }
}


void CSXReadSNV::OutputInVal2Scrn()
{
   fprintf(stderr,"Chromosome\tOffset\tRefbase\tVarBase\tQryPos\tQryScore\n");

   for (unsigned long i=0; i<m_ulTotalRecs; i++)
   {
       fprintf(stderr,"%d\t\t%u\t%c\t%c\t%d\t%d\n",m_pSNV[i].cChromosome,
       ComputeOffset(m_pSNV[i].acOffset),m_pSNV[i].cRefBase,m_pSNV[i].cVarBase,
       m_pSNV[i].ucQryPos,m_pSNV[i].ucQltyScore);
   }

   fprintf(stderr,"Total Rec Read = %u\n",m_ulTotalRecs);
}


void CSXReadSNV::OutputInVal2File()
{
    FILE *pfOut = fopen(m_pcOutFile,"w");
    if (!pfOut) {fprintf(stderr,"Can't open output file..."); return;}

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        fprintf(pfOut,"%d\t%u\t%c\t%c\t%d\t%d\n",m_pSNV[i].cChromosome,
                ComputeOffset(m_pSNV[i].acOffset),m_pSNV[i].cRefBase,
                m_pSNV[i].cVarBase,m_pSNV[i].ucQryPos,m_pSNV[i].ucQltyScore);
    }

    fclose(pfOut);
}


void CSXReadSNV::OutputSNVStats(const char c, const char* pcFile)
{
    SNV_STATS_LIST::iterator itr;
    char acReportName[100]; char acSummRptName[100]; char acType[10];

    int nMinQP=999,nMaxQP=1;
    if (c=='a')//All
    {
        sprintf(acReportName,"SNV Quality Score Analysis");
        sprintf(acSummRptName,"SNV Summary by Base Quality Score");
        sprintf(acType,"SNV");

        for (unsigned long i=0; i<m_ulTotalRecs; i++)
        {
            itr = m_SNVStats_List.find(&m_pSNV[i]);

            if (itr==m_SNVStats_List.end()){
                m_SNVStats_List[&m_pSNV[i]]++;
                if (m_pSNV[i].ucQryPos > nMaxQP) nMaxQP = m_pSNV[i].ucQryPos;
                if (m_pSNV[i].ucQryPos < nMinQP) nMinQP = m_pSNV[i].ucQryPos;
            }else{
                itr->second++;
            }
        }
    }
    else if (c=='s')//SNP
    {
        sprintf(acReportName,"SNP Quality Score Analysis");
        sprintf(acSummRptName,"SNP Summary by Base Quality Score");
        sprintf(acType,"SNP");

        for (unsigned long i=0; i<m_ulTotalRecs; i++)
        {
            if (m_pSNV[i].cRefBase!='-' && m_pSNV[i].cVarBase!='-')
            {
                itr = m_SNVStats_List.find(&m_pSNV[i]);

                if (itr==m_SNVStats_List.end()){
                    m_SNVStats_List[&m_pSNV[i]]++;
                    if (m_pSNV[i].ucQryPos > nMaxQP) nMaxQP = m_pSNV[i].ucQryPos;
                    if (m_pSNV[i].ucQryPos < nMinQP) nMinQP = m_pSNV[i].ucQryPos;
                }else{
                    itr->second++;
                }
            }
        }        
    }
    else if (c=='i')//INS
    {
        sprintf(acReportName,"INS Quality Score Analysis");
        sprintf(acSummRptName,"INS Summary by Base Quality Score");
        sprintf(acType,"INS");

        for (unsigned long i=0; i<m_ulTotalRecs; i++)
        {
            if (m_pSNV[i].cRefBase =='-')
            {
                itr = m_SNVStats_List.find(&m_pSNV[i]);

                if (itr==m_SNVStats_List.end()){
                    m_SNVStats_List[&m_pSNV[i]]++;
                    if (m_pSNV[i].ucQryPos > nMaxQP) nMaxQP = m_pSNV[i].ucQryPos;
                    if (m_pSNV[i].ucQryPos < nMinQP) nMinQP = m_pSNV[i].ucQryPos;
                }else{
                    itr->second++;
                }
            }
        }
    }
    else if (c=='d')//DEL
    {
        sprintf(acReportName,"DEL Quality Score Analysis");
        sprintf(acSummRptName,"DEL Summary by Base Quality Score");
        sprintf(acType,"DEL");

        for (unsigned long i=0; i<m_ulTotalRecs; i++)
        {
            if (m_pSNV[i].cVarBase =='-')
            {
                itr = m_SNVStats_List.find(&m_pSNV[i]);

                if (itr==m_SNVStats_List.end()){
                    m_SNVStats_List[&m_pSNV[i]]++;
                    if (m_pSNV[i].ucQryPos > nMaxQP) nMaxQP = m_pSNV[i].ucQryPos;
                    if (m_pSNV[i].ucQryPos < nMinQP) nMinQP = m_pSNV[i].ucQryPos;
                }else{
                    itr->second++;
                }
            }
        }
    }

    itr = m_SNVStats_List.begin();

    if (itr == m_SNVStats_List.end()) return;
    int nStart =(*(itr->first)).ucQltyScore;
    itr = m_SNVStats_List.end(); itr--;
    int nEnd =(*(itr->first)).ucQltyScore;

    int nTotalRows = nMaxQP-nMinQP+1;
    int nTotalCols = nEnd-nStart+1;
    int anStatsTbl[nTotalRows][nTotalCols];
    
    for (int i=0; i <nTotalRows; i++)
    {
        for (int j=0;j<nTotalCols;j++)        
           anStatsTbl[i][j] = 0;                                           
    }
    
    itr = m_SNVStats_List.begin();

    unsigned long ulAllSNPs,ulOutRangeSNPs,ulQSLess20,ulQS20,
    ulQS21,ulQS22,ulQS23,ulQS24,ulQS25,ulQS26,ulQSGreater26;

    ulAllSNPs=ulOutRangeSNPs=ulQSLess20=ulQS20=ulQS21=ulQS22=
    ulQS23=ulQS24=ulQS25=ulQS26=ulQSGreater26=0;

    int nQScore,nCol,nRow;
    while(itr!= m_SNVStats_List.end())
    {
        nRow =  (*itr->first).ucQryPos - nMinQP;
        nCol = (*itr->first).ucQltyScore - nStart;
        
        anStatsTbl[nRow][nCol] = itr->second;

        if  ((*itr->first).ucQryPos ==1||(*itr->first).ucQryPos >31) {
            ulOutRangeSNPs+=itr->second;
        }
        else{
            nQScore = (*itr->first).ucQltyScore -m_nQScoreOfs;

            if  (nQScore < 20) {ulQSLess20+=itr->second;}
            else if (nQScore == 20) {ulQS20+=itr->second;}
            else if (nQScore == 21) {ulQS21+=itr->second;}
            else if (nQScore == 22) {ulQS22+=itr->second;}
            else if (nQScore == 23) {ulQS23+=itr->second;}
            else if (nQScore == 24) {ulQS24+=itr->second;}
            else if (nQScore == 25) {ulQS25+=itr->second;}
            else if (nQScore == 26) {ulQS26+=itr->second;}
            else if (nQScore > 26) {ulQSGreater26+=itr->second;}
        }
        ulAllSNPs+=itr->second;

        itr++;
    }

    m_SNVStats_List.clear();

    FILE *pfStats=fopen(pcFile,"w");
    if (!pfStats) {fprintf(stderr,"Failed to open %s ...",pcFile); return;}

    fprintf(pfStats,"#REPORT NAME\t%s\n",acReportName);
    fprintf(pfStats,"#PROJECT NAME\n");
    fprintf(pfStats,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfStats,"#LANE NO\tALL\n");
    fprintf(pfStats,"#GENERATED AT\t"); PrintRptDateTime(pfStats);
    fprintf(pfStats,"#PROGRAM & BUILD\tReadSNV Rev %s\n", VERSION);
    fprintf(pfStats,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfStats,"#FILTER\t>=1\n\n\n");

    fprintf(pfStats,"QryPos/QScore");
    for (int j=0;j<nTotalCols;j++)
       fprintf(pfStats,"\t%d",j+nStart-m_nQScoreOfs);

    fprintf(pfStats,"\n");


    for (int i=0;i<nTotalRows;i++)
    {
        fprintf(pfStats,"%d",nMinQP++);
        for (int j=0;j<nTotalCols;j++){
            fprintf(pfStats,"\t%d",anStatsTbl[i][j]);
        }
        fprintf(pfStats,"\n");
    }

    fclose(pfStats);

    //output SNV Summary table

     struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acSummFName[strlen(m_pcSampleID)+100];
    sprintf(acSummFName,"%s_%c%c%c_stat_bqs_%02d%02d%02d.rpt",
            m_pcSampleID,tolower(acType[0]),tolower(acType[1]),tolower(acType[2]),
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    FILE *pf = fopen(acSummFName,"w");

    if (!pf){fprintf(stderr,"Failed to open %s ...",acSummFName); goto Exit;}

    unsigned long ulTotalSNPs=ulAllSNPs-ulOutRangeSNPs;
    unsigned long ulFilteredSNPs=ulQS23+ulQS24+ulQS25+ulQS26+ulQSGreater26;

    fprintf(pf,"#REPORT NAME\t%s\n",acSummRptName);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pf,"#FILTER\t>=1\n\n\n");


    fprintf(pf,"\tCriteria\tabsolute\tpercentage\n");
    fprintf(pf,"All %s\t\t%u\t100.00\n",acType,ulAllSNPs);
    fprintf(pf,"Out-of-range %ss\tqOffset=1 or ",acType);
    fprintf(pf,"qOfset>31\t%u\t%.2f\n",ulOutRangeSNPs,(float)ulOutRangeSNPs/(float)ulAllSNPs*100);
    //fprintf(pf,"Total %ss\t\t%u\t%.2f\n",acType,ulTotalSNPs,(float)ulTotalSNPs/(float)ulAllSNPs*100);
    fprintf(pf,"\n");
    //fprintf(pf,"SNP-to-read ratio(%)\t\txxx\tyyy\n");
    //fprintf(pf,"\n");
    fprintf(pf,"Base Quality Score Analysis:\n");
    fprintf(pf,"Base Quality Score (BQS)\tQS<20\t%u\t%.2f\n",ulQSLess20,(float)ulQSLess20/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=20\t%u\t%.2f\n",ulQS20,(float)ulQS20/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=21\t%u\t%.2f\n",ulQS21,(float)ulQS21/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=22\t%u\t%.2f\n",ulQS22,(float)ulQS22/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=23\t%u\t%.2f\n",ulQS23,(float)ulQS23/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=24\t%u\t%.2f\n",ulQS24,(float)ulQS24/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=25\t%u\t%.2f\n",ulQS25,(float)ulQS25/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=26\t%u\t%.2f\n",ulQS26,(float)ulQS26/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS>26\t%u\t%.2f\n",ulQSGreater26,(float)ulQSGreater26/(float)ulTotalSNPs*100);
    fprintf(pf,"Total filtered %ss (BQS >=23)\t\t%u\t%.2f\n",acType,ulFilteredSNPs,(float)ulFilteredSNPs/(float)ulTotalSNPs*100);
    //fprintf(pf,"\n");
    //fprintf(pf,"Unique SNPs with QS>=23 and\n");
    //fprintf(pf,"minimum one supporting read\t\txxx\n");
    
 Exit:
    if (pf) fclose(pf); 
}


void CSXReadSNV::OutputSNVFreq()
{
    FILE *pf=fopen(m_pcFreqFile,"w");
    if (!pf) {fprintf(stderr,"Failed to open %s ...",m_pcFreqFile); return;}

    SNV_FREQ_LIST::iterator itr;

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        stSNVFreq *pSNVFreq = new stSNVFreq(&m_pSNV[i]);

        itr = m_SNVFreq_List.find(pSNVFreq);

        if (itr==m_SNVFreq_List.end()){
            m_SNVFreq_List[pSNVFreq]++;
        }else{
            delete pSNVFreq;
            itr->second++;
        }
    }
    
    itr = m_SNVFreq_List.begin();
    
    fprintf(pf,"Chromosome\tOffset\tStrand\tRefBase\tVarBase\tFreq\n");
    while(itr!= m_SNVFreq_List.end())
    {        
        fprintf(pf,"%d\t%u\t%c\t%c\t%c\t%u\n",abs((*(itr->first)).cChromosome),
                (*(itr->first)).unOffset,(*(itr->first)).cStrand,
                (*(itr->first)).cRefBase,(*(itr->first)).cVarBase,
                itr->second);
        itr++;
    }

    itr = m_SNVFreq_List.begin();

    while (itr!= m_SNVFreq_List.end()){
        delete itr->first; itr++;
    }

    m_SNVFreq_List.clear(); fclose(pf);
}


void CSXReadSNV::OutputChroFreq()
{
    FILE *pf=fopen(m_pcChroFreqFile,"w");
    if (!pf) {fprintf(stderr,"Failed to open %s ...",m_pcChroFreqFile); return;}

    CHRO_FREQ_LIST::iterator itr;

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        itr = m_ChroFreq_List.find(&m_pSNV[i]);

        if (itr==m_ChroFreq_List.end())
            m_ChroFreq_List[&m_pSNV[i]]++;
        else
            itr->second++;
    }

    itr = m_ChroFreq_List.begin();

    fprintf(pf, "Chromosome\tFreq\n");

    while(itr!=m_ChroFreq_List.end()){
        fprintf(pf,"%d\t%u\n",abs((*itr->first).cChromosome),itr->second); itr++;
    }

    m_ChroFreq_List.clear(); fclose(pf);
}


void CSXReadSNV::OutputChroStrand()
{
    FILE *pf=fopen(m_pcChroStrandFile,"w");
    if (!pf) {fprintf(stderr,"Failed to open %s ...",m_pcChroStrandFile); return;}

    CHRO_STRAND_LIST::iterator itr;

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        stChroStrand *pChroStrand = new stChroStrand(&m_pSNV[i]);
        itr = m_ChroStrand_List.find(pChroStrand);

        if (itr==m_ChroStrand_List.end())
            m_ChroStrand_List[pChroStrand]++;
        else
        {
            delete pChroStrand;
            itr->second++;
        }
    }
        
    itr = m_ChroStrand_List.begin();

    fprintf(pf, "Chromosome\tStrand\tFreq\n");

    while(itr!=m_ChroStrand_List.end()){
        fprintf(pf,"%d\t%c\t%u\n",(*itr->first).cChromosome,(*itr->first).cStrand,
                itr->second); itr++;
    }

    itr = m_ChroStrand_List.begin();

    while (itr!= m_ChroStrand_List.end()){
        delete itr->first; itr++;
    }

    m_ChroStrand_List.clear(); fclose(pf);
}


void CSXReadSNV::OutputChromOffset()
{
    CHRO_OFFSET_FREQ_LIST m_CO_SNP_List;

    FILE *pf=fopen(m_pcChromOffsetFile,"w");
    if (!pf) {fprintf(stderr,"Failed to open %s ...",m_pcChromOffsetFile); return;}

    CHRO_OFFSET_FREQ_LIST::iterator itr;

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        stChromOffset *pChromOffset = new stChromOffset(&m_pSNV[i]);

        if (m_pSNV[i].cRefBase!='-' && m_pSNV[i].cVarBase!='-')
        {            
            itr = m_CO_SNP_List.find(pChromOffset);

            if (itr == m_CO_SNP_List.end())
                m_CO_SNP_List[pChromOffset]++;
            else
            {
                delete pChromOffset;
                itr->second++;
            }
        }
    }

    itr = m_CO_SNP_List.begin();

    fprintf(pf, "Chromosome\tOffset\tFreq\n");

    while(itr!=m_CO_SNP_List.end()){
        fprintf(pf,"%d\t%u\t%u\n",(*itr->first).cChromosome,(*itr->first).unOffset,
                itr->second); itr++;
    }

    fclose(pf); 

    itr = m_CO_SNP_List.begin();

    while (itr!= m_CO_SNP_List.end()){
        delete itr->first; itr++;
    }
    m_CO_SNP_List.clear();
}


void CSXReadSNV::OutputDensity()
{
    FILE *pfIn = fopen(m_pcInFile,"r");
    if (!pfIn) {fprintf(stderr,"Failed to open %s ...",m_pcInFile); return;}

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acbuf[5196]; char *pChr; int nChromosome;

    char acFile[100]; unsigned int unOffset=0; 
    unsigned long ulOffset=0, ulOffsetTmp; unsigned short nCnt=0;
    FILE *pf; FILE **apf = new FILE*[24];
    char acType[4]; 

    if (m_eSNV == eSNP) strcpy(acType,"snp");   
    else if (m_eSNV == eINS) strcpy(acType,"ins");
    else if (m_eSNV == eDEL) strcpy(acType,"del");
    else {fprintf(stdout,"Invalid type specified ...\n"); goto Exit;}   

    for (int j=0;j<24; j++)
    {
        if (j < 22) sprintf(acFile,"%s.%s_density.%d.bin",m_pcSampleID,acType,j+1);        
        else if (j==22) sprintf(acFile,"%s.%s_density.x.bin",m_pcSampleID,acType);
        else sprintf(acFile,"%s.%s_density.y.bin",m_pcSampleID,acType);        

        apf[j] = fopen(acFile,"wb");

        if (!apf[j]){fprintf(stderr,"Failed to open %s ...",acFile); goto Exit;}
    }    

    while(!feof(pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),pfIn)) break;
        if (acbuf[0]=='#'||acbuf[0]=='>'||acbuf[0]=='C') continue;

        pChr = strchr(acbuf,'\r');  if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n');  if (pChr) *pChr=0;

        if (acbuf[0] == 0) continue;

        pChr = strtok(acbuf,"\t");  //Chromosome
           
        if (!isalpha(pChr[0])) nChromosome = atoi(pChr);
        else if (pChr[0]=='X') nChromosome=23;
        else if (pChr[0]=='Y') nChromosome=24;              
        else nChromosome=0;

        if (nChromosome <1 || nChromosome > 24)
           {fprintf(stdout,"Invalid chromosome number %d ...",nChromosome); continue;}

        //pChr = strtok(NULL,"\t");   //GiNumber
        pChr = strtok(NULL,"\t");   //Offset

        if (pChr){
            unOffset = atol(pChr); unOffset--;

            pf = apf[nChromosome-1];
            ulOffsetTmp = unOffset/1000;

            if (ulOffset == ulOffsetTmp){
                if (nCnt < 255) nCnt++;
            }
            else {
                fseek(pf,ulOffset,SEEK_SET);
                fwrite(&nCnt,1,1,pf);
                ulOffset=ulOffsetTmp; nCnt=1;        
            }
        }
    }

Exit:
    for (int j=0;j<24; j++){
        if (apf[j]) fclose(apf[j]);
    }

    delete [] apf;
    fclose(pfIn);
}



void CSXReadSNV::OutputREADDensity1K()
{
    typedef struct _stRead
    {
        char acRead[2];
    } stRead;

    char acSrcFile[strlen(m_pcDensFilePath)+16];
    int fd; struct stat st; size_t file_size; unsigned long ulTotalRecs;
    stRead *pstRead=NULL; unsigned int unReads; int rc;

    char acDestFile[100]; FILE *pfDest=NULL;
    unsigned long ulOffset=0, ulOffsetTmp; unsigned short nCnt=0;

    for (int i=0; i<24;i++)
    {
        sprintf(acSrcFile,"%s/Chr%02d_read_den",m_pcDensFilePath,i+1);

        if (i<22) sprintf(acDestFile,"%s.read_density.%d.bin",m_pcSampleID,i+1);
        else if (i==22) sprintf(acDestFile,"%s.read_density.x.bin",m_pcSampleID);
        else sprintf(acDestFile,"%s.read_density.y.bin",m_pcSampleID);

        fd =open(acSrcFile, O_RDONLY, S_IRUSR);
        if (fd == -1) {fprintf(stderr, "Failed to open %s ...", acSrcFile); break;}

        fstat(fd, &st); file_size = (size_t)st.st_size; ulTotalRecs = file_size/2;

        pstRead=(stRead*)mmap(NULL,file_size,PROT_READ,MAP_PRIVATE,fd,0);

        if (pstRead == (void *)-1) {
           fprintf(stderr, "Error: mmap failed: %s\n", strerror(errno)); break;
        }

        pfDest = fopen(acDestFile,"w");
        for (unsigned long ulIdx=0; ulIdx < ulTotalRecs; ulIdx++){
            unReads = pstRead[ulIdx].acRead[1]; unReads<<=8;
            unReads|= pstRead[ulIdx].acRead[0];

            ulOffsetTmp = ulIdx/1000;

            if (ulOffset == ulOffsetTmp){
                if (nCnt < 65535) nCnt+=unReads;
            }
            else {
                fseek(pfDest,ulOffset*2,SEEK_SET);
                fwrite(&nCnt,2,1,pfDest);
                ulOffset=ulOffsetTmp; nCnt=unReads;
            }
       }
       fclose(pfDest);
       rc = munmap((stRead*)pstRead, file_size);
       if (rc == -1) {fprintf(stderr, "munmap %s failed",acSrcFile); break;}
   }
}


inline double CSXReadSNV::CalcAvgQS(unsigned long ulTotalQS, unsigned long ulTotalCnt)
{
    if (ulTotalCnt == 0) return 0;
    return (ulTotalQS == 0)? 0:(double)ulTotalQS/ulTotalCnt;
}


void CSXReadSNV::OutputSNVList()
{    
    FILE *pfSrc=NULL;

    pfSrc = fopen(m_pcInFile,"r"); if (!pfSrc) {fprintf(stdout,"Failed to open %s ...\n",m_pcInFile); goto Exit;}
    fprintf(stdout,"Loading Annotation Data Source...\n");

    m_pSXAnnotate = new SXAnnotate;

    if (!m_pSXAnnotate->LoadDataSource(m_pcBinnedFile,m_pcExceptionFile,m_pcAnnotateFile)) goto Exit;
    
    m_punGINums = strcmp(m_pcGenomeVer,"36.3")==0?g_unGINums_G36:g_unGINums_G37;  

    if (m_eSNV==eSNP){
        fprintf(stdout,"Printing SNP List Table...\n");PrintSNPList(pfSrc);
    }
    else if (m_eSNV==eINS){
        fprintf(stdout,"Printing INS List Table...\n");PrintINSList(pfSrc);
    }
    else if (m_eSNV==eDEL){
        fprintf(stdout,"Printing DEL List Table...\n");PrintDELList(pfSrc);
    }
    
Exit:

    if (pfSrc) fclose(pfSrc); 
    if (m_pSXAnnotate) {delete m_pSXAnnotate; m_pSXAnnotate=NULL;}
}


void CSXReadSNV::OutputType(FILE *pf,char cRefBase, char cVarBase)
{
    bool bTs=false;

    if ((cRefBase == 'A' && cVarBase =='G')|| (cRefBase == 'G' && cVarBase =='A')||
        (cRefBase == 'C' && cVarBase =='T')|| (cRefBase == 'T' && cVarBase =='C'))
    {
        fprintf(pf,"Ts"); bTs=true;
    }
    else if (cRefBase == 'A')
    {
        if (cVarBase =='C'|| cVarBase == 'T')
            fprintf(pf,"Tv");
    }
    else if (cRefBase == 'C')
    {
        if (cVarBase =='A'|| cVarBase == 'G')
           fprintf(pf,"Tv");
    }
    else if (cRefBase == 'G')
    {
        if (cVarBase =='C'|| cVarBase == 'T')
           fprintf(pf,"Tv");
    }
    else if (cRefBase =='T')
    {
        if (cVarBase =='A'|| cVarBase =='G')
        {
            fprintf(pf,"Tv");
        }
    }
    else{
        fprintf(pf,"-"); return;
    }

    if (bTs)
        m_stSNVAnalysis.ulTotalSNP_Transitions++;
    else
        m_stSNVAnalysis.ulTotalSNP_Transversions++;
}


inline void CSXReadSNV::ExtractAnnotation(char *pRecs, stAnnotate &aAnnotate)
{
    if (!pRecs[0]) return;
    char *pChr, *pcTmp = strchr(pRecs,'\n');
        
    while(pcTmp)
    {
        pChr = pRecs; pRecs=pcTmp+1; *pcTmp ='\0';

        pChr = strtok(pChr,"\t"); //chromosome
        pChr = strtok(NULL,"\t"); //Start Pos
        pChr = strtok(NULL,"\t"); //Stop Pos
        pChr = strtok(NULL,"\t"); //Feature type

        if (pChr[0]=='G' && pChr[1]=='N'){
            pChr = strtok(NULL,"\t"); //FeatureID
            pChr = strtok(NULL,"\t"); //Data Source
            pChr = strtok(NULL,"\t"); //ID/Name

            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.GeneName_List.push_back(pcTmp);

            pChr = strtok(NULL,"\t"); //PharmacogenomicID(Antibody)
            pChr = strtok(NULL,"\t"); //GeneType

            pChr = strtok(NULL,"\t"); //Description
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.GeneDesc_List.push_back(pcTmp);

            pChr = strtok(NULL,"\t"); //Keyword
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.GeneKW_List.push_back(pcTmp);
        }
        else if (pChr[0]=='E' && pChr[1]=='X'){
            pChr = strtok(NULL,"\t"); //FeatureID

            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.ExonID_List.push_back(pcTmp);
        }
        else if (pChr[0]=='C' && pChr[1]=='N'){
            pChr = strtok(NULL,"\t"); //FeatureID
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.CNV_List.push_back(pcTmp);
        }
        else if (pChr[0]=='M' && pChr[1]=='R'){
            pChr = strtok(NULL,"\t"); //FeatureID
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.miRNA_List.push_back(pcTmp);
        }
        else if (pChr[0]=='P' && pChr[1]=='R'){
            pChr = strtok(NULL,"\t"); //FeatureID
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.Promoter_List.push_back(pcTmp);
        }
        else if (pChr[0]=='U' && pChr[1]=='T'){
            pChr = strtok(NULL,"\t"); //FeatureID
            pcTmp = new char[strlen(pChr)+1]; strcpy(pcTmp,pChr);
            aAnnotate.UTR_List.push_back(pcTmp);
        }


        pcTmp = strchr(pRecs,'\n'); 
    }
}


inline void CSXReadSNV::PrintVECListContents(VEC_LIST &vec, FILE *pf)
{
    if (vec.size() > 0){
        for (unsigned int nCnt = 0; nCnt < vec.size(); nCnt++){
            if (nCnt > 0) fprintf(pf,",%s",vec[nCnt]);
            else fprintf(pf,"%s",vec[nCnt]);
        }
    }
    else{fprintf(pf,"-");}
}


void CSXReadSNV::PrintSNPList(FILE* pfSrc)
{
    char acbuf[4096],acAllele[4],acdbSnp[1024],acHomoHet[10],cChromosome=0;
    char *pChr=NULL; float fLocalCNV=0.0,fAvgQScore=0.0,fAvgQryPos=0.0;
    unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulFwdSR=0,ulRvsSR=0,ulSNPReads=0,ulReadDensity=0,ulPEnds=0;
    unsigned long ulMappableBases=0, ulMappableReadDens=0, ulMappableReptDens=0;   
    unsigned long ulMappableReptCnt=0, ulMappableReptVal=0,ulTest1=0,ulTest2=0;     
    bool bGene,bExon; int nSNPIdx=0; stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    
    FILE *pfSNP,*pfSyno,*pfSNPStats; pfSNP=pfSyno=pfSNPStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)];

    int nTotalTabs=0; char *pRecs=NULL; stAnnotate aAnnotate;

    sprintf(acFileName,"%s_snp_c1_%02d%02d%02d",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_snp_summary_c1_%02d%02d%02d.rpt",m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s.syn",m_pcSampleID);   //pfSyno = fopen("SNP.syn","w");

    pfSyno = fopen(acFileName,"w");
    if (!pfSyno) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}
        
    fprintf(pfSNP,"#REPORT NAME\tSNP\n");
    fprintf(pfSNP,"#PROJECT NAME\n");
    fprintf(pfSNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfSNP,"#LANE NO\tALL\n");
    fprintf(pfSNP,"#GENERATED AT\t"); PrintRptDateTime(pfSNP);
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"Chromosome\tGiNumber\tOffset\tNucleotide_Variant\tFwd_SNP_Reads"
                  "\tRvs_SNP_Reads\tTotal_SNP_Reads\tTotal_Read_Density\tPEnd_Count"
                  "\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases"
                  "\tMappable_Read_Densities\tMappable_Repeat_Densities"  
                  "\tMappable_Bases_Repeat\tRepeat_Densities_MBP" 
                  "\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts)"
                  "\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant"
                  "\tMissense\tNonsense\n");

    while (!feof(pfSrc))
    {                        
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;            
            fprintf(stdout,"Total Recs Processed = %u X 0.1M\n",ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' || acbuf[0]==0 || acbuf[0]=='C'){
           continue;
        }

        if (nTotalTabs == 0)
        {
           pChr = strchr(acbuf,'\t');
           while(pChr){
               nTotalTabs++; pChr = strchr(pChr+1,'\t');
           }
        } 
         
        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfSNP,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if(pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}      
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}
       
        if (nTotalTabs == 32){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulSNPReads=atol(pChr); //Total_SNP_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);//Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); // Read Strength
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr); //Total_Fwd_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr); //Total_Rvs_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr); //Total_PEnd
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr); //Total_AVQryPos     
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local_CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr); //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densies
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappale_Repeat_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr); //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr); //Repeat_Densities_MBP
        }
        else if(nTotalTabs == 12){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr);    //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr);    //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);       //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);       //Rvs_SNP_Reads 
            pChr = strtok(NULL,"\t"); ulSNPReads=atol(pChr);    //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);       //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);   //Avg_QryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);    //Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);     //dbSnp
//printf("a\n");fprintf(stdout,"acdbSnp => %s\n",pChr);

//            pChr = strtok(NULL,"\t"); //Reserved for Local CNV, this field has empty value at the moment
//printf("b\n");
            pChr = strtok(NULL,"\t"); 
//fprintf(stdout,"acHomoHet => %s\n",pChr);
            strcpy(acHomoHet,pChr); //Homo_Het
//printf("c\n"); return;
            //acHomoHet[0]=0;
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }  

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfSNP,"\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulSNPReads,ulReadDensity);
        
        fprintf(pfSNP,"%u\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

        pRecs = m_pSXAnnotate->GetRecord(cChromosome,ulOffset,ulOffset);

        ExtractAnnotation(pRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfSNP);fprintf(pfSNP,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfSNP);fprintf(pfSNP,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfSNP); fprintf(pfSNP,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfSNP); fprintf(pfSNP,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfSNP,"Y\t");
        else fprintf(pfSNP,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfSNP,"Y\t");
        else fprintf(pfSNP,"N\t");

        //PrintVECListContents(aAnnotate.UTR_List, pfSNP); fprintf(pfSNP,"\t");  

        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfSNP);fprintf(pfSNP,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfSNP,"-\t");}

        //if (acHomoHet[0]!=0) {
        if (nTotalTabs == 32){
           fprintf(pfSNP,"%s\t",acHomoHet); 
           fprintf(pfSNP,"%.2f\t%u\t%u\t%u\t%u\t%u\t",
                   fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                   ulMappableReptCnt,ulMappableReptVal);
        }
        else fprintf(pfSNP,"%s\t-\t-\t-\t-\t-\t-\t",acHomoHet); //"-\t-\t");//,acHomoHet);               
        
        if (aAnnotate.CNV_List.size() > 0) fprintf(pfSNP,"Y\t");
        else fprintf(pfSNP,"N\t");

        OutputType(pfSNP,acAllele[0],acAllele[2]);fprintf(pfSNP,"\n");

        nSNPIdx = (ulSNPReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulSNPReads-1;

        astSNPStats[nSNPIdx].unCnt++;
        astSNPStats[nSNPIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astSNPStats[nSNPIdx].undbSNP++;
        if (bGene) astSNPStats[nSNPIdx].unGene++;
        if (bExon){
            astSNPStats[nSNPIdx].unExon++;
            fprintf(pfSyno,"%u|%u|%u|%c%c|",astSNPStats[nSNPIdx].unCnt,unGINum,
                    ulOffset,acAllele[0],acAllele[2]);
            PrintVECListContents(aAnnotate.GeneName_List, pfSyno);
            fprintf(pfSyno,"\n");
       }

       aAnnotate.ClrList();
    }
   
    fprintf(stdout,"Printing SNP Statistic Report...\n");

    if (pfSNPStats) PrintStatsRpt(pfSNPStats,astSNPStats,"SNP");

    //fprintf(stdout,"Printing SNP Analysis Report...\n");
    //if (pfSNPStats) PrintSNVAnalysis();
  

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfSyno) fclose(pfSyno);
    if (pfSNPStats) fclose(pfSNPStats);
}


void CSXReadSNV::PrintINSList(FILE* pfSrc)
{   
    char acbuf[4096],acAllele[4],acdbSnp[1024],acHomoHet[10]; char *pChr=NULL, cChromosome=0;
    float fLocalCNV=0.0,fAvgQScore=0.0,fAvgQryPos=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulINSReads=0,ulFwdSR=0,ulRvsSR=0,ulReadDensity=0,ulPEnds=0;
    unsigned long ulMappableBases=0,ulMappableReadDens=0,ulMappableReptDens=0;
    unsigned long ulMappableReptCnt=0,ulMappableReptVal=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nINSIdx; stSNVTblStats astINSStats[m_unMaxSupportingReads+1];

    FILE *pfINS,*pfINSStats; pfINS=pfINSStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)]; int nTotalTabs=0;
    char *pcRecs=NULL; stAnnotate aAnnotate;

    sprintf(acFileName,"%s_ins_c1_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_ins_summary_c1_%02d%02d%02d.rpt",m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pfINS,"#REPORT NAME\tInsertion\n");
    fprintf(pfINS,"#PROJECT NAME\n");
    fprintf(pfINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfINS,"#LANE NO\tALL\n");
    fprintf(pfINS,"#GENERATED AT\t"); PrintRptDateTime(pfINS);
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfINS,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfINS,"Chromosome\tGiNumber\tOffset\tInserted_Base\tFwd_INS_Reads\tRvs_INS_Reads"
                  "\tTotal_INS_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                  "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                  "\tKnown_CNV_Region\n");

    while (!feof(pfSrc))
    {
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u X 0.1M\n",ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' ||acbuf[0]==0|| acbuf[0]=='C'){
           continue;
        }

        if (nTotalTabs == 0)
        {
           pChr = strchr(acbuf,'\t');

           while(pChr)
           {
               nTotalTabs++; pChr = strchr(pChr+1,'\t');
           }
        }

        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfINS,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if (pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}        
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        if (nTotalTabs == 32){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr);         //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr);         //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulINSReads=atol(pChr);         //Total_SNP_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr);      //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);         //Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);          //dbSnp
            pChr = strtok(NULL,"\t");                                //Read Strength
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);            //Total_Fwd_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);            //Total_Rvs_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);            //Total_PEnd
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);        //Total_AVQryPos 
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);        //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr);        //Local CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr);    //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappable_Repeat_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr);  //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr);  //Repeat_Densities_MBP 
        }
        else if(nTotalTabs == 12){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr);    //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr);    //Inserted base
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);       //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);       //Rvs_SNP_Reads
            pChr = strtok(NULL,"\t"); ulINSReads=atol(pChr);    //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);       //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);   //Avg_QryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);    //Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);     //dbSnp
            //pChr = strtok(NULL,"\t");                           //Reserved for Local CNV, this field has empty value at the moment  
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);   //Homo_Het
            //acHomoHet[0]=0;    
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfINS,"\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulINSReads,ulReadDensity);

        fprintf(pfINS,"%u\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

        pcRecs = m_pSXAnnotate->GetRecord(cChromosome,ulOffset,ulOffset);
        ExtractAnnotation(pcRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfINS);fprintf(pfINS,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfINS);fprintf(pfINS,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfINS); fprintf(pfINS,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfINS); fprintf(pfINS,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfINS,"Y\t");
        else fprintf(pfINS,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfINS,"Y\t");
        else fprintf(pfINS,"N\t");

        //PrintVECListContents(aAnnotate.UTR_List, pfINS); fprintf(pfINS,"\t");

        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfINS);fprintf(pfINS,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfINS,"-\t");}

        if (acHomoHet[0] != 0) {
           fprintf(pfINS,"%s\t",acHomoHet); 
           fprintf(pfINS,"%.2f\t%u\t%u\t%u\t%u\t%u\t",
                   fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                   ulMappableReptCnt,ulMappableReptVal);
        }
        else fprintf(pfINS,"%s\t-\t-\t-\t-\t-\t-\t",acHomoHet); //"-\t-\t");//,acHomoHet);

        if (aAnnotate.CNV_List.size() > 0) fprintf(pfINS,"Y\n");
        else fprintf(pfINS,"N\n");
     
        nINSIdx = (ulINSReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulINSReads-1;

        astINSStats[nINSIdx].unCnt++;
        astINSStats[nINSIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astINSStats[nINSIdx].undbSNP++;
        if (bGene) astINSStats[nINSIdx].unGene++;
        if (bExon) astINSStats[nINSIdx].unExon++;
        aAnnotate.ClrList();
    }

    fprintf(stdout,"Printing INS Statistic Report...\n");

    if (pfINSStats) PrintStatsRpt(pfINSStats,astINSStats,"INS");

ExitFunc:
    if (pfINS) fclose(pfINS); if (pfINSStats) fclose(pfINSStats);
}


void CSXReadSNV::PrintDELList(FILE* pfSrc)
{   
    char acbuf[4096],acAllele[4],acdbSnp[1024],acHomoHet[10]; char *pChr=NULL, cChromosome;
    float fLocalCNV=0.0,fAvgQryPos=0.0,fAvgQScore=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulFwdSR=0,ulRvsSR=0,ulDELReads=0,ulPEnds=0,ulReadDensity=0;
    unsigned long ulMappableBases=0,ulMappableReadDens=0,ulMappableReptDens=0; 
    unsigned long ulMappableReptCnt=0,ulMappableReptVal=0,ulTest1=0,ulTest2=0;

    bool bGene,bExon; int nDELIdx; stSNVTblStats astDELStats[m_unMaxSupportingReads+1];

    FILE *pfDEL,*pfDELStats; pfDEL=pfDELStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)]; int nTotalTabs=0;
    char *pcRecs=NULL; stAnnotate aAnnotate;

    sprintf(acFileName,"%s_del_c1_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_del_summary_c1_%02d%02d%02d.rpt",m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pfDEL,"#REPORT NAME\tDeletion\n");
    fprintf(pfDEL,"#PROJECT NAME\n");
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfDEL,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfDEL,"Chromosome\tGiNumber\tOffset\tDeleted_Base\tFwd_DEL_Reads\tRvs_DEL_Reads"
                  "\tTotal_DEL_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                  "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                  "\tKnown_CNV_Region\n");

    while (!feof(pfSrc))
    {
        ulTest1++;
        if (ulTest1 == 100000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %u X 0.1M\n",ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#'||acbuf[0]==0||acbuf[0]=='C'){
           continue;
        }

        if (nTotalTabs == 0)
        {
           pChr = strchr(acbuf,'\t');

           while(pChr)
           {
               nTotalTabs++; pChr = strchr(pChr+1,'\t');
           }
        }

        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfDEL,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if (pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}        
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        if (nTotalTabs == 32){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr);         //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr);         //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulDELReads=atol(pChr);         //Total_DEL_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr);      //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);         //Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);          //dbSnp
            pChr = strtok(NULL,"\t");                                //Read Strength
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);            //Total_Fwd_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);            //Total_Rvs_SNP
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);            //Total_PEnd
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);        //Total_AVQryPos     
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);        //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr);        //Local CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr);    //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappable_Rept_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr);  //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr);  //Repeat_Densities_MBP  
        }
        else if(nTotalTabs == 12){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr);    //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr);    //Deleted_Variant
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);       //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);       //Rvs_SNP_Reads
            pChr = strtok(NULL,"\t"); ulDELReads=atol(pChr);    //DEL_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);       //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);   //Avg_QryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);    //Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);     //dbSnp
            //pChr = strtok(NULL,"\t"); //Reserved for Local CNV, this field has empty value at the moment 
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);   //Homo_Het
            //acHomoHet[0]=0;  
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfDEL,"\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulDELReads,ulReadDensity);

        fprintf(pfDEL,"%u\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);  

        pcRecs = m_pSXAnnotate->GetRecord(cChromosome,ulOffset,ulOffset);
        ExtractAnnotation(pcRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfDEL);fprintf(pfDEL,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfDEL);fprintf(pfDEL,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfDEL); fprintf(pfDEL,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfDEL); fprintf(pfDEL,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfDEL,"Y\t");
        else fprintf(pfDEL,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfDEL,"Y\t");
        else fprintf(pfDEL,"N\t"); 

        //PrintVECListContents(aAnnotate.UTR_List, pfDEL); fprintf(pfDEL,"\t");

        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfDEL);fprintf(pfDEL,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfDEL,"-\t");}

        if (acHomoHet[0]!=0) {
           fprintf(pfDEL,"%s\t",acHomoHet); 
           fprintf(pfDEL,"%.2f\t%u\t%u\t%u\t%u\t%u\t",
                   fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                   ulMappableReptCnt,ulMappableReptVal);
        }
        else fprintf(pfDEL,"%s\t-\t-\t-\t-\t-\t-\t",acHomoHet); //"-\t-\t");//,acHomoHet);

        if (aAnnotate.CNV_List.size() > 0) fprintf(pfDEL,"Y\n");
        else fprintf(pfDEL,"N\n");

        nDELIdx = (ulDELReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulDELReads-1;

        astDELStats[nDELIdx].unCnt++;
        astDELStats[nDELIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astDELStats[nDELIdx].undbSNP++;
        if (bGene) astDELStats[nDELIdx].unGene++;
        if (bExon) astDELStats[nDELIdx].unExon++;
        aAnnotate.ClrList();
    }

    fprintf(stdout,"Printing DEL Statistic Report...\n");

    if (pfDELStats) PrintStatsRpt(pfDELStats,astDELStats,"DEL");

ExitFunc:
    if (pfDEL) fclose(pfDEL); if (pfDELStats) fclose(pfDELStats);
}


void CSXReadSNV::PrintStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{
    fprintf(pf,"#REPORT NAME\t%s Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);    
    fprintf(pf,"#FILTER\t%u<=Supporting Read<=%u\n\n\n",
                m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pf,"Supporting_Reads\t%s_Counts\tAvg_Read_Density\tdbSNP_Count\t"
               "dbSNP_%%\tGene_Count\tGene_%%\tExon_Count\tExon_%%\t\t\t"
               "Accum_SupportingReads\tAccum_%s_Counts\tAccum_%s_Counts_%%\t"
               "Accum_dbSNP_Count\tAccum_dbSNP_%%\tAccum_Gene_Count\t"
               "Accum_Gene_%%\tAccum_Exon_Count\tAccum_Exon_%%\n",pcType,pcType,pcType);
    //Calc Totals
    unsigned long ulMax = m_unMaxSupportingReads+1;

    for (unsigned int i=0; i < ulMax; i++)
    {
        m_stSNVAnalysis.ulTotalSNP+=pstTblStats[i].unCnt;
        m_stSNVAnalysis.ulTotalSNP_dbsnp+=pstTblStats[i].undbSNP;
        m_stSNVAnalysis.ulTotalSNP_gene+=pstTblStats[i].unGene;
        m_stSNVAnalysis.ulTotalSNP_exon+=pstTblStats[i].unExon;
    }

    fprintf(pf,"=%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t\t\t"
               ">=%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\n",
               m_unMinSupportingReads,pstTblStats[0].unCnt,
               Division(pstTblStats[0].ulTotalReadDens,pstTblStats[0].unCnt),
               pstTblStats[0].undbSNP,
               Percentage(pstTblStats[0].undbSNP,pstTblStats[0].unCnt),
               pstTblStats[0].unGene,
               Percentage(pstTblStats[0].unGene,pstTblStats[0].unCnt),
               pstTblStats[0].unExon,
               Percentage(pstTblStats[0].unExon,pstTblStats[0].unCnt),
               m_unMinSupportingReads,m_stSNVAnalysis.ulTotalSNP,
               (m_stSNVAnalysis.ulTotalSNP >0)?100.0:0.0,
               m_stSNVAnalysis.ulTotalSNP_dbsnp,
               Percentage(m_stSNVAnalysis.ulTotalSNP_dbsnp,m_stSNVAnalysis.ulTotalSNP),
               m_stSNVAnalysis.ulTotalSNP_gene,
               Percentage(m_stSNVAnalysis.ulTotalSNP_gene,m_stSNVAnalysis.ulTotalSNP),
               m_stSNVAnalysis.ulTotalSNP_exon,
               Percentage(m_stSNVAnalysis.ulTotalSNP_exon,m_stSNVAnalysis.ulTotalSNP));

    unsigned long ulAccumSnp,ulAccumdbSnp,ulAccumGene,ulAccumExon;

    ulAccumSnp=m_stSNVAnalysis.ulTotalSNP;
    ulAccumdbSnp=m_stSNVAnalysis.ulTotalSNP_dbsnp;
    ulAccumGene=m_stSNVAnalysis.ulTotalSNP_gene;
    ulAccumExon=m_stSNVAnalysis.ulTotalSNP_exon;

    for (unsigned int i=m_unMinSupportingReads; i < m_unMaxSupportingReads; i++){
        ulAccumSnp-=pstTblStats[i-1].unCnt;
        ulAccumdbSnp-=pstTblStats[i-1].undbSNP;
        ulAccumGene-=pstTblStats[i-1].unGene;
        ulAccumExon-=pstTblStats[i-1].unExon;
        
        fprintf(pf,"=%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t\t\t"
                   ">=%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\n",
                   i+1,pstTblStats[i].unCnt,
                   Division(pstTblStats[i].ulTotalReadDens,pstTblStats[i].unCnt),
                   pstTblStats[i].undbSNP,
                   Percentage(pstTblStats[i].undbSNP,pstTblStats[i].unCnt),
                   pstTblStats[i].unGene,
                   Percentage(pstTblStats[i].unGene,pstTblStats[i].unCnt),
                   pstTblStats[i].unExon,
                   Percentage(pstTblStats[i].unExon,pstTblStats[i].unCnt),
                   i+1,ulAccumSnp,Percentage(ulAccumSnp,m_stSNVAnalysis.ulTotalSNP),
                   ulAccumdbSnp,Percentage(ulAccumdbSnp,ulAccumSnp),
                   ulAccumGene,Percentage(ulAccumGene,ulAccumSnp),
                   ulAccumExon,Percentage(ulAccumExon,ulAccumSnp));
    }

    fprintf(pf,">%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t\t\t"
               ">=%u\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\t%u\t%.2f\n",
               m_unMaxSupportingReads,pstTblStats[m_unMaxSupportingReads].unCnt,
               Division(pstTblStats[m_unMaxSupportingReads].ulTotalReadDens,
                        pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].undbSNP,
               Percentage(pstTblStats[m_unMaxSupportingReads].undbSNP,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].unGene,
               Percentage(pstTblStats[m_unMaxSupportingReads].unGene,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].unExon,
               Percentage(pstTblStats[m_unMaxSupportingReads].unExon,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               m_unMaxSupportingReads+1,pstTblStats[m_unMaxSupportingReads].unCnt,
               Percentage(pstTblStats[m_unMaxSupportingReads].unCnt,
                          m_stSNVAnalysis.ulTotalSNP),
               pstTblStats[m_unMaxSupportingReads].undbSNP,
               Percentage(pstTblStats[m_unMaxSupportingReads].undbSNP,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].unGene,
               Percentage(pstTblStats[m_unMaxSupportingReads].unGene,
                          pstTblStats[m_unMaxSupportingReads].unCnt),
               pstTblStats[m_unMaxSupportingReads].unExon,
               Percentage(pstTblStats[m_unMaxSupportingReads].unExon,
                          pstTblStats[m_unMaxSupportingReads].unCnt));

    fprintf(pf,"Total\t%u\t\t%u\t\t%u\t\t%u\n",
            m_stSNVAnalysis.ulTotalSNP,m_stSNVAnalysis.ulTotalSNP_dbsnp,
            m_stSNVAnalysis.ulTotalSNP_gene,m_stSNVAnalysis.ulTotalSNP_exon);
    fprintf(pf,"\n\n\n");
}


void CSXReadSNV::PrintSNVAnalysis()
{
    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[28+strlen(m_pcSampleID)];
    sprintf(acFileName,"%s_snv_analysis_c1_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    FILE *pf = fopen(acFileName,"w");
    if (!pf) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pf,"#REPORT NAME\tSNV Analysis\n");
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#FILTER\t%u<=Supporting Read<=%u, Read Strength>=%.2f\n\n\n",
            m_unMinSupportingReads,m_unMaxSupportingReads,m_fMinReadStrength);

    fprintf(pf,"SNP Transitions/transversions ratio\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalSNP_Transitions,m_stSNVAnalysis.ulTotalSNP_Transversions));

    fprintf(pf,"SNP total count\t%u\n",m_stSNVAnalysis.ulTotalSNP);
    fprintf(pf,"INS total count\t%u\n",m_stSNVAnalysis.ulTotalINS);
    fprintf(pf,"DEL total count\t%u\n",m_stSNVAnalysis.ulTotalDEL);

    fprintf(pf,"SNP novel rate\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalSNP-m_stSNVAnalysis.ulTotalSNP_dbsnp,
                       m_stSNVAnalysis.ulTotalSNP));

    fprintf(pf,"INS novel rate\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalINS-m_stSNVAnalysis.ulTotalINS_dbsnp,
                       m_stSNVAnalysis.ulTotalINS));

    fprintf(pf,"DEL novel rate\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalDEL-m_stSNVAnalysis.ulTotalDEL_dbsnp,
                       m_stSNVAnalysis.ulTotalDEL));

    fprintf(pf,"Insertion/deletion ratio\t%.2f\n",
            Division(m_stSNVAnalysis.ulTotalINS,m_stSNVAnalysis.ulTotalDEL));

    fprintf(pf,"Ins+del/SNP ratio\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalINS+m_stSNVAnalysis.ulTotalDEL,
                       m_stSNVAnalysis.ulTotalSNP));

    fprintf(pf,"Coding Insertion/deletions ratio\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalINS_exon,m_stSNVAnalysis.ulTotalDEL_exon));

    fprintf(pf,"Coding SNP/all SNP ratio\t%.6f\n",
            Division(m_stSNVAnalysis.ulTotalSNP_exon,m_stSNVAnalysis.ulTotalSNP));

ExitFunc:
    if (pf) fclose(pf);

}


inline void CSXReadSNV::OutputChrom(FILE *pf, unsigned int unChrom)
{
    if (unChrom < 23)
        fprintf(pf,"%d",unChrom);
    else if (unChrom == 23)
        fprintf(pf,"X");
    else
        fprintf(pf,"Y");
}


unsigned char g_acBuf[2];unsigned int g_unReads;
inline unsigned int CSXReadSNV::GetReadDensity(FILE *pf, unsigned int uioffset)
{
    if (!pf) return 0;
    uioffset--;
    fseek(pf,uioffset*2,SEEK_SET);
    fread(g_acBuf,2,1,pf);

    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|= g_acBuf[0];

    return g_unReads;
}


inline unsigned int CSXReadSNV::GetReadDensityEx(FILE *pf, unsigned int uioffset)
{
    if (!pf) return 0;

    unsigned int unRead=0, unTotalReads=0;
    uioffset--;

    for (int i=0;i<2;i++)
    {
         fseek(pf,(uioffset+i)*2,SEEK_SET);
         fread(g_acBuf,2,1,pf);

         unRead = g_acBuf[1];unRead<<=8;
         unRead|= g_acBuf[0];
         unTotalReads += unRead;
    }

/*
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
 */
    g_unReads = (unsigned int)((double)((double)(unTotalReads)/2)+0.5);

    return g_unReads;
}


inline double CSXReadSNV::Percentage(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator*100;
}


inline double CSXReadSNV::Division(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator;
}


char g_acMonth[12][4]= {"Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec"};

inline void CSXReadSNV::PrintRptDateTime(FILE *pf)
{
    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    fprintf(pf,"%02d-%s-20%02d; %02d:%02d:%02d\n",
               ptm->tm_mday,g_acMonth[ptm->tm_mon],ptm->tm_year-100,
               ptm->tm_hour, ptm->tm_min, ptm->tm_sec);
}


void CSXReadSNV::TraceSNVSupportingReads()
{
    SNV_TRACE_LIST pSNVTrace_List; SNV_TRACE_LIST::iterator itr;
    QRY_POS_LIST::iterator QryPos_itr; QLY_SCORE_LIST::iterator QlyScore_itr;
    stSNVTraceData *pSNVTraceData=NULL;    

    for (unsigned long i=0; i< m_ulTotalRecs; i++)
    {
        if (m_eSNV ==eSNP)
        {
            if (m_pSNV[i].cRefBase=='-' || m_pSNV[i].cVarBase=='-')
                continue;
        }
        else if (m_eSNV ==eINS)
        {
            if (m_pSNV[i].cRefBase!='-' || m_pSNV[i].cVarBase=='-')
                continue;
        }
        else if (m_eSNV ==eDEL)
        {
            if (m_pSNV[i].cRefBase=='-' && m_pSNV[i].cVarBase!='-')
                continue;
        }

#if _TESTING
        if (i > 8000){
            //fprintf(stderr,"Total rec = %u\n",i);
         break;
        }
#endif
        stSNVTrace *pSNVTrace = new stSNVTrace(&m_pSNV[i]);

        itr = pSNVTrace_List.find(pSNVTrace);

        if (itr== pSNVTrace_List.end()){
            pSNVTraceData = new stSNVTraceData;
            pSNVTraceData->unCnt=1;

            pSNVTraceData->QryPosList[m_pSNV[i].ucQryPos]=1;
            pSNVTraceData->QlyScoreList[m_pSNV[i].ucQltyScore]=1;

            pSNVTrace_List[pSNVTrace]=pSNVTraceData;
        }
        else{
            QryPos_itr = (*itr->second).QryPosList.find(m_pSNV[i].ucQryPos);

            if (QryPos_itr==(*itr->second).QryPosList.end())
                (*itr->second).QryPosList[m_pSNV[i].ucQryPos]=1;
            else
                QryPos_itr->second++;

            QlyScore_itr = (*itr->second).QlyScoreList.find(m_pSNV[i].ucQltyScore);

            if (QlyScore_itr==(*itr->second).QlyScoreList.end())
                (*itr->second).QlyScoreList[m_pSNV[i].ucQltyScore]=1;
            else
                QlyScore_itr->second++;

            (*itr->second).unCnt++;

            delete pSNVTrace;
        }
    }

    int nQryPosCnt=0,nQlyScoreCnt=0;

    FILE *pf = fopen("SNVTrace.lst","w");

    itr = pSNVTrace_List.begin();

    while (itr != pSNVTrace_List.end())
    {
        if ((*itr->second).unCnt == m_unMinSupportingReads)
        {
            QryPos_itr = (*itr->second).QryPosList.begin();
            QlyScore_itr = (*itr->second).QlyScoreList.begin();

            for (unsigned int i=0; i < m_unMinSupportingReads; i++)
            {
                //fprintf(pf,"CNT = %u\t",(*itr->second).unCnt);
                fprintf(pf,"%u\t%u\t%c\t%c\t",m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                        m_pSNVClustered[i].cRefBase,m_pSNVClustered[i].cVarBase);

                if (QryPos_itr!= (*itr->second).QryPosList.end()){
                    fprintf(pf, "%u\t",QryPos_itr->first);
                    nQryPosCnt=QryPos_itr->second;
                }
                else{
                    fprintf(pf, "-\t");
                }

                if (--nQryPosCnt == 0) QryPos_itr++;

                if (QlyScore_itr!= (*itr->second).QlyScoreList.end()){
                    fprintf(pf, "%d\t",QlyScore_itr->first-m_nQScoreOfs);
                    nQlyScoreCnt=QlyScore_itr->second;
                }
                else{
                    fprintf(pf, "-\t");
                }

                if (--nQlyScoreCnt == 0) QlyScore_itr++;

                fprintf(pf, "\n");
            }
        }

        itr++;
    }

    fclose(pf);
}


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


void CSXReadSNV::AllocMergeListSize()
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

 
void *SortSNVData(void *ptr)
{
    stSNVRecList *pSNVRecList = (stSNVRecList*)ptr;
    char tmpChr;
//toupper(m_pSNV[i].cRefBase)
    for (int i=0; i < pSNVRecList->nUsed; i++)
    {
        pSNVRecList->pSNVRec[i].unOffset = ComputeOffset(pSNVRecList->pSNVRec[i].acOffset); 
        tmpChr = toupper(pSNVRecList->pSNVRec[i].cRefBase);
        pSNVRecList->pSNVRec[i].cRefBase = tmpChr;  
        tmpChr = toupper(pSNVRecList->pSNVRec[i].cVarBase);
        pSNVRecList->pSNVRec[i].cVarBase = tmpChr;
    } 
    qsort(pSNVRecList->pSNVRec,pSNVRecList->nUsed,sizeof(stSNVRec), CompareSNVRecs);
    return 0; 
}


inline void CSXReadSNV::GenerateSNVTbl()
{
    /*        
    FILE *pfIn=fopen("unsort.txt","w");

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
           fprintf(pfIn,"%d\t%u\t%c\t%c\t%d\t%d\t%d\t%d\n",
                m_pSNV[i].cChromosome,ComputeOffset(m_pSNV[i].acOffset),
                toupper(m_pSNV[i].cRefBase),toupper(m_pSNV[i].cVarBase),
                m_pSNV[i].ucQryPos,m_pSNV[i].ucQltyScore,m_pSNV[i].usPEnd,
                m_pSNV[i].cStrand);
    }
    fclose(pfIn);
    */

    unsigned long ulTest1=0, ulTest2=0;

#if _TESTING
    m_ulTotalRecs=25000; //Testing
#endif

    fprintf(stdout,"Loading, converting and splitting...\n"); //Converting chromosome to absolute value and compute offset...\n");

    int nChro; stSNVRecList aSNVRecList[24]; stSNVRecList *pSNVRecList=NULL;    

    FILE *pfVRI = fopen("IndenticalVRBase.lst","w");  

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {                
        if (m_pSNV[i].cRefBase == m_pSNV[i].cVarBase) {fprintf(pfVRI,"%d\n",i); continue;}
 
        nChro = abs(m_pSNV[i].cChromosome); pSNVRecList = &aSNVRecList[nChro-1];

        if (pSNVRecList->nUsed == pSNVRecList->nTotal){
            pSNVRecList->nTotal += SNVREC_SIZE;
            pSNVRecList->pSNVRec = (stSNVRec*)realloc(pSNVRecList->pSNVRec, sizeof(stSNVRec)*pSNVRecList->nTotal); 
        }

        memcpy(pSNVRecList->pSNVRec[pSNVRecList->nUsed].acOffset,m_pSNV[i].acOffset,4);
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cRefBase = m_pSNV[i].cRefBase; 
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cVarBase = m_pSNV[i].cVarBase; 
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQryPos = m_pSNV[i].ucQryPos;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQltyScore = m_pSNV[i].ucQltyScore;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].usPEnd = m_pSNV[i].usPEnd;  
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cStrand = m_pSNV[i].cStrand; 
        pSNVRecList->nUsed++;
    }

    ClrMemoryMap(); fclose(pfVRI);
    
    fprintf(stdout,"sorting...\n");

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

    /*  
    FILE *pf=fopen("sorted.txt","w");
      
        pSNVRecList = &aSNVRecList[21];//i];
        for (int j =0; j < pSNVRecList->nUsed; j++)
        {
             fprintf(pf,
                "22\t%u\t%c\t%c\t%d\t%d\t%d\t%d\n",pSNVRecList->pSNVRec[j].unOffset,
                toupper(pSNVRecList->pSNVRec[j].cRefBase),toupper(pSNVRecList->pSNVRec[j].cVarBase),
                pSNVRecList->pSNVRec[j].ucQryPos,pSNVRecList->pSNVRec[j].ucQltyScore,
                pSNVRecList->pSNVRec[j].usPEnd,pSNVRecList->pSNVRec[j].cStrand);
        }

    */

    /*
    FILE *pf=fopen("sorted.txt","w");
    for (unsigned long i=0; i<24; i++)
    {
        pSNVRecList = &aSNVRecList[i];
        for (int j =0; j < pSNVRecList->nUsed; j++)
        {
             fprintf(pf,
                "%d\t%u\t%c\t%c\t%d\t%d\t%d\t%d\n",i+1,pSNVRecList->pSNVRec[j].unOffset,
                toupper(pSNVRecList->pSNVRec[j].cRefBase),toupper(pSNVRecList->pSNVRec[j].cVarBase),
                pSNVRecList->pSNVRec[j].ucQryPos,pSNVRecList->pSNVRec[j].ucQltyScore,
                pSNVRecList->pSNVRec[j].usPEnd,pSNVRecList->pSNVRec[j].cStrand);
        }
    }
    
    fclose(pf);        
    */

    fprintf(stdout,"alloc memory for clustering...\n");

    AllocMergeListSize();

    fprintf(stdout,"clustering...\n");

    FILE *pfLog = fopen("InvalidRec.lst","w");
    char cChromosome_Curr,cRefBase_Curr,cVarBase_Curr;
    char cRef_ACGT,cVar_ACGT,cStrand;
    unsigned int unOffset_Curr,unOffset;
    unsigned short usQP, usQS; unsigned short usPEnd=0;

    int nChroIdx;
    for (nChroIdx=0; nChroIdx<24; nChroIdx++){
         if (aSNVRecList[nChroIdx].nUsed < 1) continue;
         pSNVRecList=&aSNVRecList[nChroIdx]; break; 
    }
     
    unsigned short ausPEType[4]={49152,32768,16384,0};

    for (int idx=0; idx < 4; idx++)
    {
       if ((pSNVRecList->pSNVRec[nChroIdx].usPEnd & ausPEType[idx])==ausPEType[idx]){
          usPEnd = pSNVRecList->pSNVRec[nChroIdx].usPEnd - ausPEType[idx];
          break;
       }    
    }

    //usPEnd = pSNVRecList->pSNVRec[0].usPEnd; 

    usQP = pSNVRecList->pSNVRec[0].ucQryPos;
    usQS = pSNVRecList->pSNVRec[0].ucQltyScore-m_nQScoreOfs;
    cChromosome_Curr = nChroIdx+1;      
    cStrand = pSNVRecList->pSNVRec[0].cStrand;
    unOffset_Curr = pSNVRecList->pSNVRec[0].unOffset;
    cRefBase_Curr = pSNVRecList->pSNVRec[0].cRefBase;
    cVarBase_Curr = pSNVRecList->pSNVRec[0].cVarBase;

    m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
    m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
    m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
    m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
    m_pSNVClustered[merge_list_used].ulTotalQP = usQP;  
    m_pSNVClustered[merge_list_used].ulTotalQS = usQS;
    m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
    m_pSNVClustered[merge_list_used].unCnt = 1;
    m_pSNVClustered[merge_list_used].bKeep = true;

    if (usPEnd != 0) m_pSNVClustered[merge_list_used].unPEndCnt = 1;

    if (cStrand==0) m_pSNVClustered[merge_list_used].unRvsCnt = 1;
    else if (cStrand==1) m_pSNVClustered[merge_list_used].unFwdCnt = 1;
              
    if ((cRefBase_Curr!='A' && cRefBase_Curr!='C' && cRefBase_Curr!='G' && cRefBase_Curr!='T' && 
         cRefBase_Curr!='-')|| (cVarBase_Curr!='A' && cVarBase_Curr!='C' && cVarBase_Curr!='G' &&
         cVarBase_Curr!='T' && cVarBase_Curr!='-'))
    {
        m_pSNVClustered[merge_list_used].bKeep=false;
        //LogInvalidRec(pfLog,pTmp[0]);
    }
   
   for (int i=0; i < 24; i++)
   {
       pSNVRecList=&aSNVRecList[i];
        
       for (int j=0; j < pSNVRecList->nUsed; j++)
       {   

           if (i==nChroIdx && j==0) continue;

           if ((pSNVRecList->pSNVRec[j].cRefBase!='A' && pSNVRecList->pSNVRec[j].cRefBase!='C' 
                && pSNVRecList->pSNVRec[j].cRefBase!='G'&& pSNVRecList->pSNVRec[j].cRefBase!='T' 
                && pSNVRecList->pSNVRec[j].cRefBase!='-')||
               (pSNVRecList->pSNVRec[j].cVarBase!='A' && pSNVRecList->pSNVRec[j].cVarBase!='C' 
                && pSNVRecList->pSNVRec[j].cVarBase!='G' && pSNVRecList->pSNVRec[j].cVarBase!='T' 
                && pSNVRecList->pSNVRec[j].cVarBase!='-')             
              )
           {
               //LogInvalidRec(pfLog,pTmp[i]); 
              continue;
           }
   
           cRef_ACGT = pSNVRecList->pSNVRec[j].cRefBase;
           cVar_ACGT = pSNVRecList->pSNVRec[j].cVarBase;
           unOffset = pSNVRecList->pSNVRec[j].unOffset;
           usQP = pSNVRecList->pSNVRec[j].ucQryPos;
           usQS = pSNVRecList->pSNVRec[j].ucQltyScore-m_nQScoreOfs;

           for (int idx=0; idx < 4; idx++)
    	   {
       	       if ((pSNVRecList->pSNVRec[j].usPEnd & ausPEType[idx])==ausPEType[idx])
               {
          	  usPEnd = pSNVRecList->pSNVRec[j].usPEnd - ausPEType[idx];
          	  break;
       	       }
    	   }
  
           //usPEnd = pSNVRecList->pSNVRec[j].usPEnd;

           cStrand = pSNVRecList->pSNVRec[j].cStrand;

           if (cChromosome_Curr != i+1||unOffset_Curr != unOffset||
               cRefBase_Curr != cRef_ACGT||cVarBase_Curr != cVar_ACGT) 
           {
               cChromosome_Curr = i+1; unOffset_Curr = unOffset;
               cRefBase_Curr = cRef_ACGT; cVarBase_Curr = cVar_ACGT;

               merge_list_used++;

               if (merge_list_used == merge_list_alloc) AllocMergeListSize();

               m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
               m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
               m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
               m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
               m_pSNVClustered[merge_list_used].ulTotalQP = usQP;
               m_pSNVClustered[merge_list_used].ulTotalQS = usQS;              
               m_pSNVClustered[merge_list_used].ucMaxQS = usQS;               
               m_pSNVClustered[merge_list_used].unCnt = 1;
               m_pSNVClustered[merge_list_used].bKeep = true;

               if (usPEnd != 0) m_pSNVClustered[merge_list_used].unPEndCnt = 1;

               if (cStrand==0) m_pSNVClustered[merge_list_used].unRvsCnt = 1;
               else if (cStrand==1) m_pSNVClustered[merge_list_used].unFwdCnt = 1;                
          }
          else
          {
               m_pSNVClustered[merge_list_used].ulTotalQP += usQP; 
               m_pSNVClustered[merge_list_used].ulTotalQS += usQS;

               if (m_pSNVClustered[merge_list_used].ucMaxQS < usQS)
                   m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
                 
               //if (m_pSNVClustered[merge_list_used].sMaxPEnd < usPEnd)
               //    m_pSNVClustered[merge_list_used].sMaxPEnd = usPEnd;

               m_pSNVClustered[merge_list_used].unCnt++;

               if (usPEnd != 0) m_pSNVClustered[merge_list_used].unPEndCnt++;    

               if (cStrand==0) m_pSNVClustered[merge_list_used].unRvsCnt++;
               else if (cStrand==1) m_pSNVClustered[merge_list_used].unFwdCnt++;   
          }

          ulTest1++;
          if (ulTest1 == 1000000)
          {
              ulTest2++;
              fprintf(stdout,"Total Recs gone through clustering = %uM\n", ulTest2);
              ulTest1=0;
          }
       }
    }

    if (pfLog) fclose(pfLog);
    fprintf(stdout,"Total Clustered Recs = %u\n", merge_list_used+1);    
/*        
    FILE *pfCluster=fopen("clustered.txt","w");

    for (unsigned long i=0; i < merge_list_used+1; i++)
    {
        fprintf(pfCluster,"%d\t%u\t%c>%c\t%u\t%u\t%u\t%u\n",
                m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                m_pSNVClustered[i].cRefBase,m_pSNVClustered[i].cVarBase,
                m_pSNVClustered[i].unCnt,m_pSNVClustered[i].ulTotalQS,
                m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].ucMaxQS);
    }

    fclose(pfCluster);  
*/
    return;
      
}


void CSXReadSNV::LogInvalidRec(FILE *pf,stSNVRec &Rec,int nChromosome)
{
    fprintf(pf,"%d\t%u\t%c\t%c\t%d\t%d\n",nChromosome, Rec.unOffset,
            Rec.cRefBase, Rec.cVarBase,Rec.ucQryPos,Rec.ucQltyScore);
}

void CSXReadSNV::FilterSNVTbl()
{
    
    for (unsigned long i=0; i < merge_list_used; i++)
    {
        if (m_pSNVClustered[i].ucMaxQS < m_nMinQScore )  //Ilumina - 30, CG-0
            m_pSNVClustered[i].bKeep=false;

        //if (m_pSNVClustered[i].sMaxPEnd < 1)
        //    m_pSNVClustered[i].bKeep=false;
    }
    
    
    //FILE *pf1 = fopen("duplet.txt","w");
    //FILE *pf2 = fopen("triplet.txt","w");

    unsigned long j,k;  

    for (unsigned long i=0; i < merge_list_used; i++)
    {        
        if (!m_pSNVClustered[i].bKeep) continue;

        j=i+1;

        if (!m_pSNVClustered[j].bKeep) continue;

        if (m_pSNVClustered[j].unOffset - m_pSNVClustered[i].unOffset > 1)
            continue;

        if (m_pSNVClustered[i].cChromosome!=m_pSNVClustered[j].cChromosome)
            continue;

        k=i+2;

        if (m_pSNVClustered[k].unOffset - m_pSNVClustered[j].unOffset > 1)
        {
            FilterSNVRec(m_pSNVClustered[i], m_pSNVClustered[j]); //comparing 2 set of data
            
            //PrintBaseRepeat(pf1,m_pSNVClustered[i]);
            //PrintBaseRepeat(pf1,m_pSNVClustered[j]);fprintf(pf1,"\n");

        }
        else
        {
            FilterSNVRec(m_pSNVClustered[i], m_pSNVClustered[j],
                         m_pSNVClustered[k]);                     //comaparing 3 set of data

            //PrintBaseRepeat(pf2,m_pSNVClustered[i]);
            //PrintBaseRepeat(pf2,m_pSNVClustered[j]);
            //PrintBaseRepeat(pf2,m_pSNVClustered[k]);fprintf(pf2,"\n");
        }
    }

    //if (pf1) fclose(pf1); if (pf2) fclose(pf2);
}


void CSXReadSNV::FilterSNVRec(stSNVClusterRec &Rec1, stSNVClusterRec &Rec2)
{
    if (Rec1.unCnt < Rec2.unCnt)
    {
       Rec1.bKeep = false;
    }
    else if (Rec2.unCnt < Rec1.unCnt)
    {
       Rec2.bKeep = false;
    }
    else if  (Rec1.cRefBase != '-' && Rec1.cVarBase != '-') //SNP - low priority, don't keep
    {
        if (!(Rec2.cRefBase != '-' && Rec2.cVarBase != '-'))  
            Rec1.bKeep = false;
        else  //filter the same type
        {
            Rec1.bKeep = false;Rec2.bKeep = false;
        }
    }
    else if  (Rec2.cRefBase != '-' && Rec2.cVarBase != '-') //SNP - low priority, don't keep
    {
        if (!(Rec1.cRefBase != '-' && Rec1.cVarBase != '-'))  
            Rec2.bKeep = false;
    }
    else if  (Rec1.cVarBase == '-') //DEL - high priority, keep, thus throw the other one
    {
        if (Rec2.cVarBase != '-')  
            Rec2.bKeep = false;
        else  //filter the same type
        { 
            Rec1.bKeep = false;Rec2.bKeep = false;  
        }
    }
    else if  (Rec2.cVarBase == '-') //DEL - high priority, keep, thus throw the other one
    {
        if (Rec1.cVarBase != '-')  //do not filter the same type
            Rec1.bKeep = false;
    }
    else if  (Rec1.cRefBase == '-' && Rec2.cRefBase == '-') //Both INS - filter the same type
    {
         Rec1.bKeep = false;Rec2.bKeep = false;
    }
}


void CSXReadSNV::FilterSNVRec(stSNVClusterRec &Rec1, stSNVClusterRec &Rec2,
                              stSNVClusterRec &Rec3)
{
    if (Rec1.unCnt > Rec2.unCnt && Rec1.unCnt > Rec3.unCnt)
    {
        Rec2.bKeep = false;  Rec3.bKeep = false;
    }
    else if (Rec2.unCnt > Rec1.unCnt && Rec2.unCnt > Rec3.unCnt)
    {
        Rec1.bKeep = false; Rec3.bKeep = false;
    }
    else if (Rec3.unCnt > Rec1.unCnt && Rec3.unCnt > Rec2.unCnt)
    {
        Rec1.bKeep = false; Rec2.bKeep = false;
    }
    else if (Rec1.unCnt < Rec2.unCnt && Rec1.unCnt < Rec3.unCnt)
    {
        Rec1.bKeep = false;  FilterSNVRec(Rec2, Rec3);
    }
    else if (Rec2.unCnt < Rec1.unCnt && Rec2.unCnt < Rec3.unCnt)
    {
        Rec2.bKeep = false; FilterSNVRec(Rec1, Rec3);
    }
    else if (Rec3.unCnt < Rec1.unCnt && Rec3.unCnt < Rec2.unCnt)
    {
        Rec3.bKeep = false; FilterSNVRec(Rec1, Rec2);
    }
    else // All have same supporting reads
    {
        // All are of the same type
        if ((Rec1.cRefBase == '-' && Rec2.cRefBase == '-' && Rec3.cRefBase == '-')||
            (Rec1.cVarBase == '-' && Rec2.cVarBase == '-' && Rec3.cVarBase == '-')||
            ((Rec1.cRefBase != '-' && Rec1.cVarBase != '-') &&
             (Rec2.cRefBase != '-' && Rec2.cVarBase != '-') &&
             (Rec3.cRefBase != '-' && Rec3.cVarBase != '-')))
        {

            Rec1.bKeep = false;Rec2.bKeep = false;Rec3.bKeep = false;
            return;
        }
        else
        {
            FilterSNVRec(Rec1, Rec2);
            if (!Rec1.bKeep && Rec2.bKeep)
                FilterSNVRec(Rec2, Rec3);
            else if ((Rec1.bKeep && !Rec2.bKeep))
                FilterSNVRec(Rec1, Rec3);
            else
            {
                FilterSNVRec(Rec1, Rec3);
                if (!Rec1.bKeep)
                    Rec2.bKeep = false;
            }
        }
    }
}


void CSXReadSNV::PrintBaseRepeat(FILE *pf,stSNVClusterRec &Rec)
{

    fprintf(pf,"%d\t%u\t%c\t%c\t%u\t%s\n",Rec.cChromosome, Rec.unOffset,
            Rec.cRefBase, Rec.cVarBase,Rec.unCnt, Rec.bKeep? "Keep":"Throw");

}


void CSXReadSNV::CloseHsRefFiles()
{
    for (short i=0;i<24;i++) {
        if (m_apfHsRef[i]){ fclose(m_apfHsRef[i]); m_apfHsRef[i]=NULL;}
    }
}


void CSXReadSNV::LoadHRefSeq()
{
    struct stat st; size_t file_size;
    char acbuf[130]; char *pChr=NULL;
    int nLen,nTotalRecs; unsigned int unIdx;
    
    for (short i=0;i<24;i++)
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
       fprintf(stdout,"Human Ref. File %d => ",i+1);
       for (int j=0; j < 70; j++)
       {
            fprintf(stdout,"%c",m_apcHsRef[i][j]);
       }
       fprintf(stdout,"...\n");
//#endif
   }   
   CloseHsRefFiles();
}


void CSXReadSNV::OutputAvgQScoreTable()
{
    FILE *pfdbSnp=NULL, *pfdbIndel=NULL; char acFileName[1024];

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

    m_apfHsRef = new FILE*[24];
    m_apcHsRef = new char*[24];

    for (short i=0;i<24;i++)
    {
        if (i < 22){
            sprintf(acFileName,"%s/hs_ref_chr%d.fa",m_pcHsRefFilePath,i+1);
         }
        else if(i == 22)
            sprintf(acFileName,"%s/hs_ref_chrX.fa",m_pcHsRefFilePath);
        else
            sprintf(acFileName,"%s/hs_ref_chrY.fa",m_pcHsRefFilePath);

        m_apfHsRef[i] = fopen(acFileName,"r");

        if (m_apfHsRef[i]){
           fprintf(stdout,"Touch %s ...\n",acFileName);
        }
        else{
            fprintf(stdout,"Failed to open %s ....\n",acFileName); goto Exit;            
        }
    }

    //fprintf(stdout,"Generating dbSNP Map...\n");
    //if (!GenerateDBSNPMap(pfdbSnp)) goto Exit;
    fprintf(stdout,"Generating dbINDEL Map...\n");
    if (!GenerateDBINDELMap(pfdbIndel)) goto Exit;
    fprintf(stdout,"Generating SNV Table...\n"); GenerateSNVTbl();
    fprintf(stdout,"Filtering SNV Table...\n"); FilterSNVTbl();
    fprintf(stdout,"Load Human Reference Sequence ...\n"); LoadHRefSeq();
    fprintf(stdout,"Printing SNP, INS, DEL Lists...\n"); PrintAvgQScoreList();
   
Exit:

    if (pfdbSnp) fclose(pfdbSnp); if (pfdbIndel) fclose(pfdbIndel);
    ClrDBSNPMap(); ClrDBINDELMap();

    for (int j=0;j<24; j++){
        if (m_apfRD[j]) fclose(m_apfRD[j]);
    }

    delete[] m_apfRD;
   
    CloseHsRefFiles();
}


int CSXReadSNV::CheckPoly(int nChroIdx, unsigned int unOffset, char cAllele,eSNVType eSNV,FILE *pf)
{
    fprintf(pf,"%d\t%d\t=>\t",nChroIdx+1,unOffset);

    if (eSNV == eDEL) unOffset--; // if it is an INS we need to shift 1 bp to the right. Thus, no need to minus the offset by 1
       
    if (m_apcHsRef[nChroIdx][unOffset]!=cAllele) { 
       if (eSNV == eDEL)unOffset++; 
       fprintf(pf,"%d\t%c\n",unOffset,cAllele);
       return unOffset;
    }

    unsigned int i;   

    for (i=unOffset+1;i < m_anHsRefSize[nChroIdx]; i++)
    {
        if (m_apcHsRef[nChroIdx][i] != cAllele)
           break; 
    }

/////Testing///////////////////////////////
    if (unOffset != i-1){
       int nSize = i - (unOffset+1); 
       char acAllele[nSize+1];
       memcpy(&acAllele[0],&m_apcHsRef[nChroIdx][unOffset],nSize);
       acAllele[nSize]=0;
       fprintf(pf,"%d\t%c\t%s\n",i,cAllele,acAllele);      
    }
    else
    {
       fprintf(pf,"*%d\t%c\t%c\n",i,cAllele,m_apcHsRef[nChroIdx][unOffset]);
    }  

///////////////////////////////////////////

   return i;
}


void CSXReadSNV::AddReadCnt(unsigned long ulIdxTo, unsigned long ulFrom)
{
    m_pSNVClustered[ulIdxTo].unCnt += m_pSNVClustered[ulFrom].unCnt;
    m_pSNVClustered[ulIdxTo].unFwdCnt += m_pSNVClustered[ulFrom].unFwdCnt;
    m_pSNVClustered[ulIdxTo].unRvsCnt += m_pSNVClustered[ulFrom].unRvsCnt;
    m_pSNVClustered[ulIdxTo].ulTotalQS += m_pSNVClustered[ulFrom].ulTotalQS;
    m_pSNVClustered[ulIdxTo].ulTotalQP += m_pSNVClustered[ulFrom].ulTotalQP;
}


void CSXReadSNV::PrintAvgQScoreList()
{
/*
    FILE *pfPoly = fopen("Poly.lst","w");
    FILE *pfPolySNP = fopen("PolySNP.lst","w");
*/
    char acPolyINSFileName[1024]; sprintf(acPolyINSFileName,"PolyINS_%s.lst",m_pcSampleID);
    FILE *pfPolyINS = fopen(acPolyINSFileName,"w");

    char acPolyDELFileName[1024]; sprintf(acPolyDELFileName,"PolyDEL_%s.lst",m_pcSampleID);
    FILE *pfPolyDEL = fopen(acPolyDELFileName,"w"); 

    char acThrowINSFileName[1024]; sprintf(acThrowINSFileName,"ThrowINS_%s.lst",m_pcSampleID);
    FILE *pfThrowINS = fopen(acThrowINSFileName,"w");

    char acThrowDELFileName[1024]; sprintf(acThrowDELFileName,"ThrowDEL_%s.lst",m_pcSampleID);
    FILE *pfThrowDEL = fopen(acThrowDELFileName,"w");

    char acRDFileName[1024]; sprintf(acRDFileName,"PolyRD_DEL_%s.lst",m_pcSampleID);
    FILE *pfPolyRD_DEL = fopen(acRDFileName,"w");

    //unsigned int unTmpRD=0;
/////////////////////////////

    stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astDELStats[m_unMaxSupportingReads+1];
    int nSNPIdx,nINSIdx,nDELIdx=0;   

    unsigned long ulTest1=0, ulTest2=0;
    stChkPoly ChkPoly_INS, ChkPoly_DEL;
    bool /*bSNP1=true,*/ bINS1=true, bDEL1=true;
    unsigned int unChroIdx=0,unReadDensity=0,unOffset=0;

    bool bdbSNP=false; char acRsid[200];

    float fAvgQryPos,fAvgQScore;

    FILE *pfSNP,*pfINS,*pfDEL,*pfSyno,*pfSNPStats,*pfINSStats,*pfDELStats;
    pfSNP=pfINS=pfDEL=pfSyno=pfSNPStats=0,pfINSStats=pfDELStats=0;

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[24+strlen(m_pcSampleID)];
    sprintf(acFileName,"%s_snp_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}


    sprintf(acFileName,"%s_snp_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    sprintf(acFileName,"%s_ins_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_ins_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////
    sprintf(acFileName,"%s_del_c1_aqs_%02d%02d%02d.lst",m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s_del_summary_c1_aqs_%02d%02d%02d.rpt",m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    fprintf(pfSNP,"#REPORT NAME\tSNP Avg. Score Distribution\n");
    fprintf(pfSNP,"#PROJECT NAME\n");
    fprintf(pfSNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfSNP,"#LANE NO\tALL\n");
    fprintf(pfSNP,"#GENERATED AT\t"); PrintRptDateTime(pfSNP);
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfSNP,"#FILTER\t>=1\n\n\n");
    //fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n\n\n",
    //              m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"Chromosome\tOffset\tNucleotide_Variant\tFwd_SNP_Reads"
                  "\tRvs_SNP_Reads\tTotal_SNP_Reads" 
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");

    fprintf(pfINS,"#REPORT NAME\tINS Avg. Score Distribution\n");
    fprintf(pfINS,"#PROJECT NAME\n");
    fprintf(pfINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfINS,"#LANE NO\tALL\n");
    fprintf(pfINS,"#GENERATED AT\t"); PrintRptDateTime(pfINS);
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfINS,"#FILTER\t>=1\n\n\n");

    fprintf(pfINS,"Chromosome\tOffset\tInserted_Base\tFwd_INS_Reads"
                  "\tRvs_INS_Reads\tTotal_INS_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");

    fprintf(pfDEL,"#REPORT NAME\tDEL Avg. Score Distribution\n");
    fprintf(pfDEL,"#PROJECT NAME\n");
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfDEL,"Chromosome\tOffset\tDeleted_Base\tFwd_DEL_Reads"                  
                  "\tRvs_DEL_Reads\tTotal_DEL_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");
 
    fprintf(pfPolyRD_DEL,"Chromosome\tNew_Offset1\tOrg_Offset1\tNew_RD1\tOrg_RD1\t"
                         "Org_SR\tNew_RD2\tOrg_RD2\n");

    for (unsigned long i=0; i < (merge_list_used+1); i++)
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Clustered Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }
        
        if (m_pSNVClustered[i].cChromosome == 0) continue;

        unChroIdx = m_pSNVClustered[i].cChromosome-1;

        if (m_pSNVClustered[i].cRefBase!='-' && m_pSNVClustered[i].cVarBase!='-')
        {
            bdbSNP=OutputSNPNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                  m_pSNVClustered[i].cVarBase,pfSNP,acRsid);
             
            //if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);            
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue; 
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfSNP,m_pSNVClustered[i].cChromosome);
            fprintf(pfSNP,"\t%u\t%c>%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,
                    m_pSNVClustered[i].unCnt);

            fprintf(pfSNP,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfSNP,"\n");

            nSNPIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astSNPStats[nSNPIdx].unCnt++;
            if (bdbSNP) astSNPStats[nSNPIdx].undbSNP++;
            astSNPStats[nSNPIdx].ulTotalReadDens += unReadDensity;
          
        }
        else if (m_pSNVClustered[i].cRefBase =='-' )
        {
            if (bINS1){
                bINS1=false;ChkPoly_INS.ulIdx=i;

                //fprintf(pfPolyINS,">>ChkPoly_INS.ulIdx=%u,%u\ti=%u,%u\n",
                //       ChkPoly_INS.ulIdx,m_pSNVClustered[ChkPoly_INS.ulIdx].unOffset,i,m_pSNVClustered[i].unOffset);

                ChkPoly_INS.unOffset=CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,eINS,pfPolyINS);
                continue;
            }

            //fprintf(pfPolyINS,">>ChkPoly_INS.ulIdx=%u,%u\ti=%u,%u\n",
            //       ChkPoly_INS.ulIdx,ChkPoly_INS.unOffset,i,m_pSNVClustered[i].unOffset);

            unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,eINS,pfPolyINS);            

            //fprintf(pfPolyINS,">>unOffset=%u\n",unOffset);

            if (m_pSNVClustered[i].cChromosome==m_pSNVClustered[ChkPoly_INS.ulIdx].cChromosome && 
                unOffset == ChkPoly_INS.unOffset && 
                m_pSNVClustered[i].cVarBase==m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase)// && (ChkPoly_INS.ulIdx !=i))
            { 
               if (ChkPoly_INS.ulIdx !=i){ AddReadCnt(i,ChkPoly_INS.ulIdx); ChkPoly_INS.ulIdx=i;}
               continue; 
            }
                
            bdbSNP=OutputINDELNovel(m_pSNVClustered[ChkPoly_INS.ulIdx].cChromosome,ChkPoly_INS.unOffset,
                                    m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase, pfINS,acRsid);                                              
            /* 
            if (!m_pSNVClustered[ChkPoly_INS.ulIdx].bKeep && !bdbSNP) 
            { 
                fprintf(pfThrowINS,"%u\t%u\t%c\n",m_pSNVClustered[ChkPoly_INS.ulIdx].cChromosome,
                                   ChkPoly_INS.unOffset,m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase); 
                goto EndINS;
            } 
            */                       

            if (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP ) goto EndINS;

            unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx],ChkPoly_INS.unOffset);

            if ((Percentage(m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt,unReadDensity)<m_fMinReadStrength)&& !bdbSNP) goto EndINS;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

            OutputChrom(pfINS,m_pSNVClustered[ChkPoly_INS.ulIdx].cChromosome);

            fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",ChkPoly_INS.unOffset,m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase,
                    m_pSNVClustered[ChkPoly_INS.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_INS.ulIdx].unRvsCnt,
                    m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

            fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                    unReadDensity,m_pSNVClustered[ChkPoly_INS.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);
            
            fprintf(pfINS,"\n");

            nINSIdx = (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt-1;

            astINSStats[nINSIdx].unCnt++; if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
            astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;
EndINS:
            ChkPoly_INS.ulIdx = i; ChkPoly_INS.unOffset = unOffset;  
        }
        else
        {
	    if (bDEL1){
                bDEL1=false;ChkPoly_DEL.ulIdx=i;
                ChkPoly_DEL.unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,eDEL,pfPolyDEL);
                ChkPoly_DEL.unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset); 

                if (ChkPoly_DEL.unOffset != m_pSNVClustered[i].unOffset) 
                {
                    //
                    //unTmpRD = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset);
                    //OutputChrom(pfPolyRD_DEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome);
                    //  
                    //fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t%u\t",ChkPoly_DEL.unOffset,
                    //                     m_pSNVClustered[i].unOffset,ChkPoly_DEL.unReadDensity,unTmpRD);
                    // 
                    if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt) ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;
                    else ChkPoly_DEL.unReadDensity = 0;
                    // 
                    //fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t1\n",m_pSNVClustered[i].unCnt,
                    //                     ChkPoly_DEL.unReadDensity,(unTmpRD+m_pSNVClustered[i].unCnt));
                    //
                } 

                continue;
            }

            //if (i < ChkPoly_DEL.ulIdx) continue;

            unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,eDEL,pfPolyDEL);    
            
            if (m_pSNVClustered[i].cChromosome==m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome && 
                unOffset == ChkPoly_DEL.unOffset &&
                m_pSNVClustered[i].cRefBase==m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase)
            {
               if (ChkPoly_DEL.ulIdx!=i) {AddReadCnt(i,ChkPoly_DEL.ulIdx);ChkPoly_DEL.ulIdx=i;}           
 //
 //              unTmpRD = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset);
 //              OutputChrom(pfPolyRD_DEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome);
 //
 //              fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t%u\t",m_pSNVClustered[ChkPoly_DEL.ulIdx].unOffset,
 //                                   m_pSNVClustered[i].unOffset,ChkPoly_DEL.unReadDensity,unTmpRD);
 //
               if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt) ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;
               else ChkPoly_DEL.unReadDensity = 0;
//
//               fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t2\n",m_pSNVClustered[i].unCnt,
//                                    ChkPoly_DEL.unReadDensity,(unTmpRD+m_pSNVClustered[i].unCnt));
//
               continue;
            }

            bdbSNP=OutputINDELNovel(m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome,ChkPoly_DEL.unOffset,
                                    m_pSNVClustered[ChkPoly_DEL.ulIdx].cVarBase, pfDEL,acRsid); 
 
           // if (!m_pSNVClustered[ChkPoly_DEL.ulIdx].bKeep && !bdbSNP) 
           //    {fprintf(pfThrowDEL,"%u\t%u\n",m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome,ChkPoly_DEL.unOffset); goto EndDEL;}

            if (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP) goto EndDEL;

	    unReadDensity = ChkPoly_DEL.unReadDensity + m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt;

            if (Percentage(m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) goto EndDEL;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

            OutputChrom(pfDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome);
            fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",ChkPoly_DEL.unOffset,m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase,
                    m_pSNVClustered[ChkPoly_DEL.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_DEL.ulIdx].unRvsCnt,
                    m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
            
            fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                    unReadDensity,m_pSNVClustered[ChkPoly_DEL.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfDEL,"\n");

            nDELIdx = (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt-1;

            astDELStats[nDELIdx].unCnt++; if (bdbSNP) astDELStats[nDELIdx].undbSNP++; 
            astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;
EndDEL:
            ChkPoly_DEL.ulIdx = i;
            ChkPoly_DEL.unOffset = unOffset;
            ChkPoly_DEL.unReadDensity = GetReadDensity(m_apfRD[unChroIdx],m_pSNVClustered[i].unOffset);
             
            if (unOffset != m_pSNVClustered[i].unOffset)
            {
               //  
               //unTmpRD = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset);
               //OutputChrom(pfPolyRD_DEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome);

               //fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t%u\t",ChkPoly_DEL.unOffset,
               //                     m_pSNVClustered[i].unOffset,ChkPoly_DEL.unReadDensity,unTmpRD);               
               //
               if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt) ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;
               else ChkPoly_DEL.unReadDensity = 0;
               // 
               //fprintf(pfPolyRD_DEL,"\t%u\t%u\t%u\t3\n",m_pSNVClustered[i].unCnt,
               //                     ChkPoly_DEL.unReadDensity,(unTmpRD+m_pSNVClustered[i].unCnt));
               //
            }             
        }
    }

//////////////////////

    if (!m_pSNVClustered[ChkPoly_INS.ulIdx].bKeep && !bdbSNP) goto ProcLastDELRec;

    fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);
    fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

    unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx], ChkPoly_INS.unOffset);
 
    if ((m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads) && !bdbSNP) goto ProcLastDELRec;
    if ((Percentage(m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt,unReadDensity) < m_fMinReadStrength) &&
                    !bdbSNP) goto ProcLastDELRec;

    OutputChrom(pfINS,m_pSNVClustered[ChkPoly_INS.ulIdx].cChromosome);

    fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",
                   ChkPoly_INS.unOffset,m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase,
                   m_pSNVClustered[ChkPoly_INS.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_INS.ulIdx].unRvsCnt,
                   m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

    fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                   unReadDensity,m_pSNVClustered[ChkPoly_INS.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

    fprintf(pfINS,"\n");

    nINSIdx = (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads)?
               m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt-1;

    astINSStats[nINSIdx].unCnt++;
    if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
    astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;

ProcLastDELRec:

    if (!m_pSNVClustered[ChkPoly_DEL.ulIdx].bKeep && !bdbSNP) goto PrintStatsRpt;

    fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
    fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

    unReadDensity = ChkPoly_DEL.unReadDensity + m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt; 

    if (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP) goto PrintStatsRpt;
    if (Percentage(m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) goto PrintStatsRpt;

    OutputChrom(pfDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].cChromosome);
    fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",
                   ChkPoly_DEL.unOffset,m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase,
                   m_pSNVClustered[ChkPoly_DEL.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_DEL.ulIdx].unRvsCnt,
                   m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

    fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                   unReadDensity,m_pSNVClustered[ChkPoly_DEL.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

    fprintf(pfDEL,"\n");

    nDELIdx = (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads)?
               m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt-1;

    astDELStats[nDELIdx].unCnt++;

    if (bdbSNP) astDELStats[nDELIdx].undbSNP++;
    astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;

//////////////////////

PrintStatsRpt:

   fprintf(stderr,"Printing AvgQScore Statistic Report...\n");

    PrintAvgQSStatsRpt(pfSNPStats,astSNPStats,"SNP");
    PrintAvgQSStatsRpt(pfINSStats,astINSStats,"INS");
    PrintAvgQSStatsRpt(pfDELStats,astDELStats,"DEL");

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfINS) fclose(pfINS);
    if (pfDEL) fclose(pfDEL); if (pfSNPStats) fclose(pfSNPStats);
    if (pfINSStats) fclose(pfINSStats); if (pfDELStats) fclose(pfDELStats);
    
//////////////////////////////////////////////////////////////////////////////
/*
    if (pfPoly) fclose(pfPoly);
    if (pfPolySNP) fclose(pfPolySNP);
*/
    if (pfPolyINS) fclose(pfPolyINS);
    if (pfPolyDEL) fclose(pfPolyDEL);

    if (pfThrowINS) fclose(pfThrowINS);
    if (pfThrowDEL) fclose(pfThrowDEL);
    if (pfPolyRD_DEL) fclose(pfPolyRD_DEL);
}



void CSXReadSNV::PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{   
    fprintf(pf,"#REPORT NAME\t%s Avg. QScore Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
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








