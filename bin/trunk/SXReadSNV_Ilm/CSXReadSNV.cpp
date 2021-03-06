/* 
 * File:   CSXReadSNV.cc
 * Author: Manish
 * 
 * Created on November 17, 2009, 12:40 PM
 */

#include "CSXReadSNV.h"
#include "SXDBIndelgChecker_OdbSnpfmt.h"

#define VERSION "13.2.0.1"

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
   m_pcSNPStatsFile=m_pcINSStatsFile=m_pcDELStatsFile=
   m_pcChroStrandFile=m_pcChromOffsetFile=m_pcSNVTbl=m_pcNewFile=
   m_pcBinnedFile=m_pcAnnotateFile=m_pcdbSnp=m_pcdbINS=m_pcdbDEL=
   m_pcdbIndel=m_pcSampleID=m_pcProjSampleID=0;

   m_punGINums=0;

   m_bOut2Screen=m_bFQryPos=m_bFQltyScore=m_bMemMapped=m_bStateRpt=m_bSNVLst=m_bFilter=
   m_bDensity=m_bSNVTrace=m_bAvgQScoreLst=m_bAvgQScoreLstMNV=m_bWAmbiguityFiles=
   m_bOutClusteredFile=m_bMNVLst=false; 
   m_nFQryPos1=m_nFQryPos2=m_nFQltyScore1=m_nFQltyScore2=-999;
   m_unMinSupportingReads=0; m_unMaxSupportingReads=70; m_unMaxReadDensity = 100;
  
   m_fMinReadStrength=0;
   m_nQScoreOfs=33; 
   m_eSNV = eALL;  m_eMNV = eMALL;

   m_pSXAnnotate = NULL;
   m_pObjdbINS=NULL;
   m_pObjdbDEL=NULL;   

   INIT_MERGED_LIST_SIZE = 500000000;
   EXT_MERGED_LIST_SIZE = 100000000;
}


CSXReadSNV::~CSXReadSNV() {
    if (m_pSXAnnotate) delete m_pSXAnnotate;
    if (m_pObjdbINS) delete m_pObjdbINS;
    if (m_pObjdbDEL) delete m_pObjdbDEL;
}

void CSXReadSNV::OutputVersion()
{
    fprintf(stdout,"Version - %s\n",VERSION);
}

void CSXReadSNV::setFilterParams(int nQryPos1,int nQryPos2,int nQScoreOfs,int nQltyScore1,
                                 int nQltyScore2,char cSNVType,const char* pcFile)
{
    m_nFQryPos1 = nQryPos1-1;
    m_nFQryPos2 = nQryPos2+1;
    m_nQScoreOfs = nQScoreOfs; 
    m_nFQltyScore1 = nQltyScore1-1;
    m_nFQltyScore2 = nQltyScore2+1;
    m_pcNewFile = pcFile;

    m_eSNV = (eSNVType)cSNVType;
}



void CSXReadSNV::setQryPos(int nQryPos1, int nQryPos2, const char* pcFile)
{
    m_bFQryPos = true;
    m_nFQryPos1 = nQryPos1-1;
    m_nFQryPos2 = nQryPos2+1;
    m_pcNewFile = pcFile;
}



void CSXReadSNV::setQltyScore(int nQScoreOfs,int nQltyScore1,int nQltyScore2,const char* pcFile)
{
    m_bFQltyScore = true;
    m_nQScoreOfs = nQScoreOfs;
    m_nFQltyScore1 = nQltyScore1-1;
    m_nFQltyScore2 = nQltyScore2+1;
    m_pcNewFile = pcFile;
}


void CSXReadSNV::setAvgQSLstInputs(const char *pcdbSnp, const char *pcdbINS, const char *pcdbDEL,
                                   const char *pcDensFilePath, unsigned int unMinSupportingReads,
                                   unsigned int unMaxSupportingReads, float fMinReadStrength,
                                   int nMinQScore, const char *pcSampleID,const char *pcGenomeVer,
                                   int nQScoreOfs,const char *pcHsRefFilePath,const char *pcOutputDir,
                                   const char* pcProjSampleID)
{
    m_pcdbSnp=pcdbSnp; m_pcdbINS=pcdbINS; m_pcdbDEL=pcdbDEL; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads; m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength=fMinReadStrength; m_nMinQScore=nMinQScore; m_pcSampleID=pcSampleID; 
    m_pcGenomeVer=pcGenomeVer; m_nQScoreOfs=nQScoreOfs; m_pcHsRefFilePath=pcHsRefFilePath;
    m_pcOutputDir=pcOutputDir; m_pcProjSampleID=pcProjSampleID;
    m_bAvgQScoreLst=true;
}


void CSXReadSNV::setInputs_gLaO(const char *pcdbSnp, const char *pcdbIndel, const char *pcDensFilePath,
                                unsigned int unMinSupportingReads, unsigned int unMaxSupportingReads,
                                float fMinReadStrength, int nMinQScore, const char *pcSampleID,
                                const char *pcGenomeVer, int nQScoreOfs, const char *pcHsRefFilePath,
                                const char *pcOutputDir,const char* pcProjSampleID)
{
    m_pcdbSnp=pcdbSnp; m_pcdbIndel=pcdbIndel; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads; m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength=fMinReadStrength; m_nMinQScore=nMinQScore; m_pcSampleID=pcSampleID;
    m_pcGenomeVer=pcGenomeVer; m_nQScoreOfs=nQScoreOfs; m_pcHsRefFilePath=pcHsRefFilePath;
    m_pcOutputDir=pcOutputDir; m_pcProjSampleID=pcProjSampleID;
    m_bAvgQScoreLst=true; 
}


void CSXReadSNV::setAvgQSLstInputsMNV(const char *pcdbSnp, const char *pcdbINS, const char *pcdbDEL,
                                      const char *pcDensFilePath,
                                      unsigned int unMinSupportingReads, unsigned int unMaxSupportingReads, 
                                      float fMinReadStrength, int nMinQScore, const char *pcSampleID,
                                      const char *pcGenomeVer, int nQScoreOfs, 
                                      const char *pcOutputDir, const char* pcProjSampleID)
{
    m_pcdbSnp=pcdbSnp; m_pcdbINS=pcdbINS; m_pcdbDEL=pcdbDEL; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads; m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength=fMinReadStrength; m_nMinQScore=nMinQScore; m_pcSampleID=pcSampleID;
    m_pcGenomeVer=pcGenomeVer; m_nQScoreOfs=nQScoreOfs; 
    m_pcOutputDir=pcOutputDir; m_pcProjSampleID=pcProjSampleID;
    m_bAvgQScoreLstMNV=true;
}


void CSXReadSNV::setSNVLstInputs(const char *pcBinnedFile, const char *pcExceptionFile,
                                 const char *pcAnnotateFile,
                                 unsigned int unMinSupportingReads,
                                 unsigned int unMaxSupportingReads,
                                 char cSNVType, const char *pcSampleID,
                                 const char *pcGenomeVer,const char *pcOutputDir,
                                 const char* pcProjSampleID)
{
    m_pcBinnedFile=pcBinnedFile; m_pcExceptionFile=pcExceptionFile;
    m_pcAnnotateFile=pcAnnotateFile;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;    
    m_eSNV = (eSNVType)cSNVType;
    m_pcSampleID=pcSampleID; 
    m_bSNVLst = true;
    m_pcGenomeVer=pcGenomeVer; m_pcOutputDir=pcOutputDir;m_pcProjSampleID=pcProjSampleID;
}


void CSXReadSNV::setMNVLstInputs(const char *pcBinnedFile, const char *pcExceptionFile,
                                 const char *pcAnnotateFile,
                                 unsigned int unMinSupportingReads,
                                 unsigned int unMaxSupportingReads,
                                 char cMNVType, const char *pcSampleID,
                                 const char *pcGenomeVer,const char *pcOutputDir,
                                 const char* pcProjSampleID)
{
    m_pcBinnedFile=pcBinnedFile; m_pcExceptionFile=pcExceptionFile;
    m_pcAnnotateFile=pcAnnotateFile;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;
    m_eMNV = (eMNVType)cMNVType;
    m_pcSampleID=pcSampleID;
    m_bMNVLst = true;
    m_pcGenomeVer=pcGenomeVer; m_pcOutputDir=pcOutputDir;m_pcProjSampleID=pcProjSampleID;
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
    if (m_bMNVLst) {OutputMNVList(); goto Exit;}
    if (m_bAvgQScoreLstMNV){OutputAvgQScoreTableMNV(); goto Exit;} 

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
    if (m_bDensity){OutputDensity();} 

    if (m_pcOutFile) {fprintf(stderr,"Generating SNV ascii file ....\n");OutputInVal2File();}
    if (m_bOut2Screen) OutputInVal2Scrn();

    if (m_pcFreqFile) {fprintf(stderr,"Generating SNV Freq Report ....\n");OutputSNVFreq();}
    if (m_pcChroFreqFile) {fprintf(stderr,"Generating Chromosome Freq Report ....\n");OutputChroFreq();}
    if (m_pcChroStrandFile) {fprintf(stderr,"Generating Chromosome Strand Freq Report ....\n");OutputChroStrand();}
    if (m_pcChromOffsetFile) {fprintf(stderr,"Generating Chromosome Offset Freq Report ....\n");OutputChromOffset();}

ClrMMap:
    ClrMemoryMap();

Exit:
    if (m_pSNVClustered) free(m_pSNVClustered);
    else  if (m_pMNVClustered) free(m_pMNVClustered);
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
        if (m_pSNV[i].nChromosome == 0) continue;

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
       if (m_pSNV[i].nChromosome == 0) continue;

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
       if (m_pSNV[i].nChromosome == 0) continue;

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
        if (m_pSNV[i].nChromosome == 0)
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
        if (m_pSNV[i].nChromosome == 0) continue;

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
        if (m_pSNV[i].nChromosome == 0) continue;

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
        if (m_pSNV[i].nChromosome == 0) continue;

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
       fprintf(stderr,"%d\t\t%u\t%c\t%c\t%d\t%d\n",m_pSNV[i].nChromosome,
       ComputeOffset(m_pSNV[i].acOffset),m_pSNV[i].cRefBase,m_pSNV[i].cVarBase,
       m_pSNV[i].ucQryPos,m_pSNV[i].ucQltyScore);
   }

   fprintf(stderr,"Total Rec Read = %lu\n",m_ulTotalRecs);
}


void CSXReadSNV::OutputInVal2File()
{
    FILE *pfOut = fopen(m_pcOutFile,"w");
    if (!pfOut) {fprintf(stderr,"Can't open output file..."); return;}

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        fprintf(pfOut,"%d\t%u\t%c\t%c\t%d\t%d\n",m_pSNV[i].nChromosome,
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

    unsigned long ulTotalSNPs=ulAllSNPs-ulOutRangeSNPs;
    unsigned long ulFilteredSNPs=ulQS23+ulQS24+ulQS25+ulQS26+ulQSGreater26;

    FILE *pf = fopen(acSummFName,"w");
    if (!pf){fprintf(stderr,"Failed to open %s ...",acSummFName); goto Exit;}

    fprintf(pf,"#REPORT NAME\t%s\n",acSummRptName);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pf,"#FILTER\t>=1\n\n\n");

    fprintf(pf,"\tCriteria\tabsolute\tpercentage\n");
    fprintf(pf,"All %s\t\t%lu\t100.00\n",acType,ulAllSNPs);
    fprintf(pf,"Out-of-range %ss\tqOffset=1 or ",acType);
    fprintf(pf,"qOfset>31\t%lu\t%.2f\n",ulOutRangeSNPs,(float)ulOutRangeSNPs/(float)ulAllSNPs*100);
    fprintf(pf,"\n");
    fprintf(pf,"Base Quality Score Analysis:\n");
    fprintf(pf,"Base Quality Score (BQS)\tQS<20\t%lu\t%.2f\n",ulQSLess20,(float)ulQSLess20/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=20\t%lu\t%.2f\n",ulQS20,(float)ulQS20/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=21\t%lu\t%.2f\n",ulQS21,(float)ulQS21/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=22\t%lu\t%.2f\n",ulQS22,(float)ulQS22/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=23\t%lu\t%.2f\n",ulQS23,(float)ulQS23/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=24\t%lu\t%.2f\n",ulQS24,(float)ulQS24/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=25\t%lu\t%.2f\n",ulQS25,(float)ulQS25/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS=26\t%lu\t%.2f\n",ulQS26,(float)ulQS26/(float)ulTotalSNPs*100);
    fprintf(pf,"Base Quality Score (BQS)\tQS>26\t%lu\t%.2f\n",ulQSGreater26,(float)ulQSGreater26/(float)ulTotalSNPs*100);
    fprintf(pf,"Total filtered %ss (BQS >=23)\t\t%lu\t%.2f\n",acType,ulFilteredSNPs,(float)ulFilteredSNPs/(float)ulTotalSNPs*100);
    
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
        fprintf(pf,"%d\t%u\t%c\t%c\t%c\t%lu\n",(*(itr->first)).nChromosome,
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
        fprintf(pf,"%d\t%lu\n",(*itr->first).nChromosome,itr->second); itr++;
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
        fprintf(pf,"%d\t%c\t%lu\n",(*itr->first).nChromosome,(*itr->first).cStrand,
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
        fprintf(pf,"%d\t%u\t%lu\n",(*itr->first).nChromosome,(*itr->first).unOffset,
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
    
    m_punGINums = strcmp(m_pcGenomeVer,"37.1")==0?g_unGINums_G37:g_unGINums_G36;  

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

    sprintf(acFileName,"%s/%s_snp_c1_%02d%02d%02d",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_snp_summary_c1_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s.syn",m_pcOutputDir,m_pcSampleID);   //pfSyno = fopen("SNP.syn","w");

    pfSyno = fopen(acFileName,"w");
    if (!pfSyno) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}
        
    fprintf(pfSNP,"#REPORT NAME\tSNP\n");
    fprintf(pfSNP,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfSNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfSNP,"#LANE NO\tALL\n");
    fprintf(pfSNP,"#GENERATED AT\t"); PrintRptDateTime(pfSNP);
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    //fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n",
    //              m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"\n\n\n");

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
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;} 
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' || acbuf[0]==0 || acbuf[0]=='C') continue;

        if (nTotalTabs == 0){ 
           pChr = strchr(acbuf,'\t'); 
           while(pChr){ nTotalTabs++; pChr = strchr(pChr+1,'\t');}
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
        else if(nTotalTabs == 12){ //||nTotalTabs == 24){
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

            pChr = strtok(NULL,"\t"); 
            strcpy(acHomoHet,pChr); //Homo_Het
        }
        else if(nTotalTabs == 17||nTotalTabs == 24){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr); //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr); //Rvs_SNP_Reads
            pChr = strtok(NULL,"\t"); ulSNPReads=atol(pChr); //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr); //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr); //AVQryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);//Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local_CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr); //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densies
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappale_Repeat_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr); //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr); //Repeat_Densities_MBP
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }  

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfSNP,"\t%u\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulSNPReads,ulReadDensity);
        
        fprintf(pfSNP,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

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

        if (acHomoHet[0] != 0){
           fprintf(pfSNP,"%s\t",acHomoHet); 
           fprintf(pfSNP,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
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
            fprintf(pfSyno,"%u|%u|%lu|%c%c|",astSNPStats[nSNPIdx].unCnt,unGINum,
                    ulOffset,acAllele[0],acAllele[2]);
            PrintVECListContents(aAnnotate.GeneName_List, pfSyno);
            fprintf(pfSyno,"\n");
       }

       aAnnotate.ClrList();

       if (strcmp(acHomoHet,"het")==0) m_stSNVAnalysis.ulTotalSNPHet++;
       else if (strcmp(acHomoHet,"hom")==0) m_stSNVAnalysis.ulTotalSNPHom++;
       else m_stSNVAnalysis.ulTotalSNPGreyArea++;
    }
   
    fprintf(stdout,"Printing SNP Statistic Report...\n");

    if (pfSNPStats) PrintStatsRpt(pfSNPStats,astSNPStats,(char*)"SNP");

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

    sprintf(acFileName,"%s/%s_ins_c1_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_ins_summary_c1_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    m_stSNVAnalysis.Init(); 

    fprintf(pfINS,"#REPORT NAME\tSingle Nucleotide Insertion\n");
    fprintf(pfINS,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfINS,"#LANE NO\tALL\n");
    fprintf(pfINS,"#GENERATED AT\t"); PrintRptDateTime(pfINS);
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    //fprintf(pfINS,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
    //              m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfINS,"\n\n\n");

    fprintf(pfINS,"Chromosome\tGiNumber\tOffset\tInserted_Base\tFwd_INS_Reads\tRvs_INS_Reads"
                  "\tTotal_INS_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                  "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                  "\tKnown_CNV_Region\n");

    while (!feof(pfSrc))
    {
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2); ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;} pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' ||acbuf[0]==0|| acbuf[0]=='C') continue;        

        if (nTotalTabs == 0) {
           pChr = strchr(acbuf,'\t');
           while(pChr){ nTotalTabs++; pChr = strchr(pChr+1,'\t');}
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
        else if(nTotalTabs == 12){ //||nTotalTabs == 24){
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
        else if(nTotalTabs == 17||nTotalTabs == 24){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr); //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr); //Rvs_SNP_Reads
            pChr = strtok(NULL,"\t"); ulINSReads=atol(pChr); //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr); //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr); //AVQryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);//Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local_CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr); //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densies
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappale_Repeat_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr); //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr); //Repeat_Densities_MBP
        }  
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfINS,"\t%u\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulINSReads,ulReadDensity);

        fprintf(pfINS,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

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
           fprintf(pfINS,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
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

        if (strcmp(acHomoHet,"het")==0) m_stSNVAnalysis.ulTotalSNPHet++;
        else if (strcmp(acHomoHet,"hom")==0) m_stSNVAnalysis.ulTotalSNPHom++;
        else m_stSNVAnalysis.ulTotalSNPGreyArea++;
    }

    fprintf(stdout,"Printing INS Statistic Report...\n");

    if (pfINSStats) PrintStatsRpt(pfINSStats,astINSStats,(char*)"INS");

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

    sprintf(acFileName,"%s/%s_del_c1_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_del_summary_c1_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pfDEL,"#REPORT NAME\tSingle Nucleotide Deletion\n");
    fprintf(pfDEL,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    //fprintf(pfDEL,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
    //              m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfDEL,"\n\n\n");

    fprintf(pfDEL,"Chromosome\tGiNumber\tOffset\tDeleted_Base\tFwd_DEL_Reads\tRvs_DEL_Reads"
                  "\tTotal_DEL_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                  "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                  "\tKnown_CNV_Region\n");

    m_stSNVAnalysis.Init();

    while (!feof(pfSrc))
    {
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#'||acbuf[0]==0||acbuf[0]=='C') continue;        

        if (nTotalTabs == 0){
           pChr = strchr(acbuf,'\t');
           while(pChr){nTotalTabs++; pChr = strchr(pChr+1,'\t');}
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
        else if(nTotalTabs == 17||nTotalTabs == 24){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr); //Fwd_SNP_Reads
            pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr); //Rvs_SNP_Reads
            pChr = strtok(NULL,"\t"); ulDELReads=atol(pChr); //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr); //PEnd_Count
            pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr); //AVQryPos
            pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);//Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local_CNV
            pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr); //Mappable_Bases
            pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr); //Mappable_Read_Densies
            pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr); //Mappale_Repeat_Densities
            pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr); //Mappable_Bases_Repeat
            pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr); //Repeat_Densities_MBP
        }   
        else if(nTotalTabs == 12){ //||nTotalTabs == 24){
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
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);   //Homo_Het
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfDEL,"\t%u\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulOffset,acAllele,ulFwdSR,ulRvsSR,ulDELReads,ulReadDensity);

        fprintf(pfDEL,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);  

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
           fprintf(pfDEL,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
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

        if (strcmp(acHomoHet,"het")==0) m_stSNVAnalysis.ulTotalSNPHet++;
        else if (strcmp(acHomoHet,"hom")==0) m_stSNVAnalysis.ulTotalSNPHom++;
        else m_stSNVAnalysis.ulTotalSNPGreyArea++;
    }

    fprintf(stdout,"Printing DEL Statistic Report...\n");

    if (pfDELStats) PrintStatsRpt(pfDELStats,astDELStats,(char*)"DEL");

ExitFunc:
    if (pfDEL) fclose(pfDEL); if (pfDELStats) fclose(pfDELStats);
}


void CSXReadSNV::PrintStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{
    fprintf(pf,"#REPORT NAME\t%s Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);    
    fprintf(pf,"#DISPLAY RANGE\t%u<=Supporting Read<=%u\n\n\n",
                m_unMinSupportingReads,m_unMaxSupportingReads);


    fprintf(pf,"Heterozygous(het)\tHomozygous(hom)\tGrey Area\tRatio het/hom\n");
    fprintf(pf,"%10lu\t%10lu\t%10lu\t%.2f\n\n\n",
            m_stSNVAnalysis.ulTotalSNPHet,m_stSNVAnalysis.ulTotalSNPHom,
            m_stSNVAnalysis.ulTotalSNPGreyArea,
            Division(m_stSNVAnalysis.ulTotalSNPHet,m_stSNVAnalysis.ulTotalSNPHom));

    if (strcmp(pcType,"SNP") == 0) {
       fprintf(pf,"Transitions(Ts)\tTransversions(Tv)\tRatio Ts/Tv\n");
       fprintf(pf,"%10lu\t%10lu\t%.2f\n\n\n",
       m_stSNVAnalysis.ulTotalSNP_Transitions,m_stSNVAnalysis.ulTotalSNP_Transversions,
       Division(m_stSNVAnalysis.ulTotalSNP_Transitions,m_stSNVAnalysis.ulTotalSNP_Transversions));
    }

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
               ">=%u\t%lu\t%.2f\t%lu\t%.2f\t%lu\t%.2f\t%lu\t%.2f\n",
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
                   ">=%u\t%lu\t%.2f\t%lu\t%.2f\t%lu\t%.2f\t%lu\t%.2f\n",
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

    fprintf(pf,"Total\t%lu\t\t%lu\t\t%lu\t\t%lu\n",
            m_stSNVAnalysis.ulTotalSNP,m_stSNVAnalysis.ulTotalSNP_dbsnp,
            m_stSNVAnalysis.ulTotalSNP_gene,m_stSNVAnalysis.ulTotalSNP_exon);
    fprintf(pf,"\n\n\n");


    if (strcmp(pcType,"SNP") != 0) return;

    fprintf(pf,"Heterozygous(het)\tHomozygous(hom)\tGrey Area\tRatio (het/hom)\n");
    fprintf(pf,"%10lu\t%10lu\t%10lu\t%.2f\n",
            m_stSNVAnalysis.ulTotalSNPHet,m_stSNVAnalysis.ulTotalSNPHom,
            m_stSNVAnalysis.ulTotalSNPGreyArea, 
            Division(m_stSNVAnalysis.ulTotalSNPHet,m_stSNVAnalysis.ulTotalSNPHom));   

    fprintf(pf,"Transitions(Ts)\tTransversions(Tv)\tRatio Ts/Tv\n");
    fprintf(pf,"%10lu\t%10lu\t%.2f\n",
            m_stSNVAnalysis.ulTotalSNP_Transitions,m_stSNVAnalysis.ulTotalSNP_Transversions,
            Division(m_stSNVAnalysis.ulTotalSNP_Transitions,m_stSNVAnalysis.ulTotalSNP_Transversions));  
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


extern "C" int CompareSNVRecs(const void *a, const void *b)
{
    if (((const stSNVRec *) a)->unOffset != ((const stSNVRec *) b)->unOffset)
        return ((const stSNVRec *) a)->unOffset - ((const stSNVRec *) b)->unOffset;

    if (((const stSNVRec *) a)->cRefBase != ((const stSNVRec *) b)->cRefBase)
        return ((const stSNVRec *) a)->cRefBase - ((const stSNVRec *) b)->cRefBase;

    return ((const stSNVRec *) a)->cVarBase - ((const stSNVRec *) b)->cVarBase;
}


void CSXReadSNV::AllocMergeListSize(stSNVClusterRec **pSNVClustered, size_t &merge_list_alloc, size_t &merge_list_used)
{
        if (*pSNVClustered == NULL) {
                merge_list_alloc = INIT_MERGED_LIST_SIZE;
                merge_list_used = 0;
                try {
                        *pSNVClustered = (stSNVClusterRec*)malloc(sizeof(stSNVClusterRec)*merge_list_alloc);
                }
                catch (...) {
                        fprintf(stderr, "Error: insufficient memory for clustered SNV reords list!\n");
                        fprintf(stderr, "Program terminated.\n");
                        exit(1);
                }

                return;
        }

        merge_list_alloc += EXT_MERGED_LIST_SIZE;
        try {
                *pSNVClustered = (stSNVClusterRec*)realloc(*pSNVClustered,sizeof(stSNVClusterRec)*merge_list_alloc);
        }
        catch (...) {
                fprintf(stderr, "Error: insufficient memory for clustered SNV reords list!\n");
                fprintf(stderr, "Program terminated.\n");
                exit(1);
        }
}


void CSXReadSNV::AllocMergeListSize()
{
	if (m_pSNVClustered == NULL) {
		m_merge_list_alloc = INIT_MERGED_LIST_SIZE;
		m_merge_list_used = 0;
		try {
                        m_pSNVClustered = (stSNVClusterRec*)malloc(sizeof(stSNVClusterRec)*m_merge_list_alloc); 
		}
		catch (...) {
			fprintf(stderr, "Error: insufficient memory for Pair Info list!\n");
			fprintf(stderr, "Program terminated.\n");
			exit(1);
		}

		return;
	}

	m_merge_list_alloc += EXT_MERGED_LIST_SIZE;
	try {
                m_pSNVClustered = (stSNVClusterRec*)realloc(m_pSNVClustered,sizeof(stSNVClusterRec)*m_merge_list_alloc);
	}
	catch (...) {
		fprintf(stderr, "Error: insufficient memory for Pair Info list\n");
		fprintf(stderr, "Program terminated.\n");
		exit(1);
	}
}

 
void *SortSNVData(void *ptr)
{
    char tmpChr; stSNVRecList *pSNVRecList = (stSNVRecList*)ptr;

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


void CSXReadSNV::PrintUnSortedSNV(stSNV *pSNV, unsigned long ulTotalRecs)
{
    fprintf(stdout,"Printing Unsorted SNV recs...\n");

    char acFileName[1024]; sprintf(acFileName,"%s/%s_unsorted.txt",m_pcOutputDir,m_pcSampleID);
    FILE *pfIn=fopen(acFileName,"w");
    unsigned long ulTest1=0, ulTest2=0;  
    
    for (unsigned long i=0; i<ulTotalRecs; i++)
    {
        if (++ulTest1==100000){fprintf(stdout,"Total Unsorted Recs being printed = %lu X 1.0M\n",++ulTest2);ulTest1=0;}
        fprintf(pfIn,"%d\t%u\t%c\t%c\t%d\t%d\t%d\t%d\n",
                pSNV[i].nChromosome,ComputeOffset(pSNV[i].acOffset),
                pSNV[i].cRefBase,pSNV[i].cVarBase,
                pSNV[i].ucQryPos,pSNV[i].ucQltyScore,pSNV[i].usPEnd,
                pSNV[i].cStrand);
    }
    fclose(pfIn);
}


void CSXReadSNV::PrintSortedSNV(stSNVRecList aSNVRecList[])
{
    fprintf(stdout,"Printing Sorted SNV recs...\n");

    char acFileName[1024]; sprintf(acFileName,"%s/%s_sorted.txt",m_pcOutputDir,m_pcSampleID);
    FILE *pf=fopen(acFileName,"w"); stSNVRecList *pSNVRecList;
    unsigned long ulTest1=0, ulTest2=0;

    for (unsigned long i=0; i<24; i++)
    {
        pSNVRecList = &aSNVRecList[i];
        for (int j =0; j < pSNVRecList->nUsed; j++)
        {
            if (++ulTest1==100000){fprintf(stdout,"Total Sorted Recs being printed = %lu X 1.0M\n",++ulTest2);ulTest1=0;}

            fprintf(pf,
                "%lu\t%u\t%c\t%c\t%d\t%d\t%d\t%d\n",i+1,pSNVRecList->pSNVRec[j].unOffset,
                toupper(pSNVRecList->pSNVRec[j].cRefBase),toupper(pSNVRecList->pSNVRec[j].cVarBase),
                pSNVRecList->pSNVRec[j].ucQryPos,pSNVRecList->pSNVRec[j].ucQltyScore,
                pSNVRecList->pSNVRec[j].usPEnd,pSNVRecList->pSNVRec[j].cStrand);
        }
    }
    
    fclose(pf);
}


void CSXReadSNV::GenerateSNVTblNew(stSNV *pSNV,unsigned long ulTotalRecs, size_t mapped_file_size,
                                   stSNVClusterRec **lpSNVClustered, unsigned long aulTotalRec[],
                                   size_t &merge_list_alloc, size_t &merge_list_used)
{
    //PrintUnSortedSNV(pSNV, ulTotalRecs); 

    fprintf(stdout,"Loading, converting and splitting...\n");
      
    unsigned long ulTest1=0, ulTest2=0;
    int nChro; stSNVRecList aSNVRecList[24]; stSNVRecList *pSNVRecList=NULL;

    for (unsigned long i=0; i<ulTotalRecs; i++)
    {
        nChro = pSNV[i].nChromosome; pSNVRecList = &aSNVRecList[nChro-1];

        if (pSNVRecList->nUsed == pSNVRecList->nTotal){
            pSNVRecList->nTotal += SNVREC_SIZE;
            pSNVRecList->pSNVRec = (stSNVRec*)realloc(pSNVRecList->pSNVRec, sizeof(stSNVRec)*pSNVRecList->nTotal);
        }

        memcpy(pSNVRecList->pSNVRec[pSNVRecList->nUsed].acOffset,pSNV[i].acOffset,4);
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cRefBase = pSNV[i].cRefBase;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cVarBase = pSNV[i].cVarBase;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQryPos = pSNV[i].ucQryPos;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].ucQltyScore = pSNV[i].ucQltyScore;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].usPEnd = pSNV[i].usPEnd;
        pSNVRecList->pSNVRec[pSNVRecList->nUsed].cStrand = pSNV[i].cStrand;
        pSNVRecList->nUsed++;
    }

    if (mapped_file_size==0) ClrMemoryMap(); //SNVList
    else{                                    //Ambiguity List
        if (m_bWAmbiguityFiles){  
            int rc = munmap((char *)pSNV,mapped_file_size);
            if (rc == -1) fprintf(stderr, "munmap failed...\n");
        }
    }

    fprintf(stdout,"sorting...\n");

    for (int i=0; i<24; i++){
         if (pthread_create(&(aSNVRecList[i].ThreadID),NULL,SortSNVData,(void*)(&aSNVRecList[i])) != 0)
         { fprintf(stdout,"Failed to spawn SortSNVData Thread %d ...\n",i); exit(0); }
    }

    for (int i=0; i<24; i++){
         if (pthread_join(aSNVRecList[i].ThreadID,NULL)!=0)
         { fprintf(stdout,"Failed to perform Thread Join %d ...\n",i); exit(0); }
         aulTotalRec[i]=0;
    }

    //PrintSortedSNV(aSNVRecList);

    fprintf(stdout,"alloc memory for clustering...\n");
    AllocMergeListSize(lpSNVClustered,merge_list_alloc,merge_list_used);

    fprintf(stdout,"clustering...\n");

    char cRefBase_Curr,cVarBase_Curr,cRef_ACGT=0,cVar_ACGT=0,cStrand;
    int nChromosome_Curr;
    unsigned int unOffset_Curr,unOffset=0;
    unsigned short usQP, usQS, usPEnd=0,ausPEType[4]={49152,32768,16384,0};

    int nChroIdx;
    for (nChroIdx=0; nChroIdx<24; nChroIdx++){
         if (aSNVRecList[nChroIdx].nUsed < 1) continue;
         pSNVRecList=&aSNVRecList[nChroIdx]; break;
    }

    for (int idx=0; idx < 4; idx++)
    {
       if ((pSNVRecList->pSNVRec[nChroIdx].usPEnd & ausPEType[idx])==ausPEType[idx])
       { usPEnd = pSNVRecList->pSNVRec[nChroIdx].usPEnd - ausPEType[idx]; break; }
    }

    usQP = pSNVRecList->pSNVRec[0].ucQryPos;

    if (pSNVRecList->pSNVRec[0].ucQltyScore < m_nQScoreOfs) usQS = 0;
    else usQS = pSNVRecList->pSNVRec[0].ucQltyScore-m_nQScoreOfs;

    nChromosome_Curr = nChroIdx+1;
    cStrand = pSNVRecList->pSNVRec[0].cStrand & 0x01;
    unOffset_Curr = pSNVRecList->pSNVRec[0].unOffset;
    cRefBase_Curr = pSNVRecList->pSNVRec[0].cRefBase;
    cVarBase_Curr = pSNVRecList->pSNVRec[0].cVarBase;

    (*lpSNVClustered)[merge_list_used].nChromosome = nChromosome_Curr;
    (*lpSNVClustered)[merge_list_used].unOffset = unOffset_Curr;
    (*lpSNVClustered)[merge_list_used].cRefBase = cRefBase_Curr;
    (*lpSNVClustered)[merge_list_used].cVarBase = cVarBase_Curr;
    (*lpSNVClustered)[merge_list_used].ulTotalQP = usQP;
    (*lpSNVClustered)[merge_list_used].ulTotalQS = usQS;
    (*lpSNVClustered)[merge_list_used].ucMaxQS = usQS;
    (*lpSNVClustered)[merge_list_used].unCnt = 1;
    (*lpSNVClustered)[merge_list_used].bKeep = true;

    aulTotalRec[0]++;

    if (usPEnd != 0) (*lpSNVClustered)[merge_list_used].unPEndCnt = 1;

    if (cStrand==0) (*lpSNVClustered)[merge_list_used].unRvsCnt = 1;
    else if (cStrand==1) (*lpSNVClustered)[merge_list_used].unFwdCnt = 1;

    if ((cRefBase_Curr!='A' && cRefBase_Curr!='C' && cRefBase_Curr!='G' && cRefBase_Curr!='T' &&
         cRefBase_Curr!='-')|| (cVarBase_Curr!='A' && cVarBase_Curr!='C' && cVarBase_Curr!='G' &&
         cVarBase_Curr!='T' && cVarBase_Curr!='-'))
    {
        (*lpSNVClustered)[merge_list_used].bKeep=false;
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
              continue;
           }

           cRef_ACGT = pSNVRecList->pSNVRec[j].cRefBase;
           cVar_ACGT = pSNVRecList->pSNVRec[j].cVarBase;
           unOffset = pSNVRecList->pSNVRec[j].unOffset;
           usQP = pSNVRecList->pSNVRec[j].ucQryPos;

           if (pSNVRecList->pSNVRec[j].ucQltyScore < m_nQScoreOfs) usQS = 0;
           else usQS = pSNVRecList->pSNVRec[j].ucQltyScore-m_nQScoreOfs;

           for (int idx=0; idx < 4; idx++)
           {
               if ((pSNVRecList->pSNVRec[j].usPEnd & ausPEType[idx])==ausPEType[idx])
               { usPEnd = pSNVRecList->pSNVRec[j].usPEnd - ausPEType[idx]; break; }
           }

	   //usPEnd = pSNVRecList->pSNVRec[j].usPEnd;
	   
           cStrand = pSNVRecList->pSNVRec[j].cStrand & 0x01;

           if (nChromosome_Curr != i+1||unOffset_Curr != unOffset||
               cRefBase_Curr != cRef_ACGT||cVarBase_Curr != cVar_ACGT)
           {
               nChromosome_Curr = i+1; unOffset_Curr = unOffset;
               cRefBase_Curr = cRef_ACGT; cVarBase_Curr = cVar_ACGT;

               merge_list_used++;

               if (merge_list_used == merge_list_alloc)
                   AllocMergeListSize(lpSNVClustered,merge_list_alloc,merge_list_used);

               (*lpSNVClustered)[merge_list_used].nChromosome = nChromosome_Curr;               
               (*lpSNVClustered)[merge_list_used].unOffset = unOffset_Curr;
               (*lpSNVClustered)[merge_list_used].cRefBase = cRefBase_Curr;
               (*lpSNVClustered)[merge_list_used].cVarBase = cVarBase_Curr;
               (*lpSNVClustered)[merge_list_used].ulTotalQP = usQP;
               (*lpSNVClustered)[merge_list_used].ulTotalQS = usQS;
               (*lpSNVClustered)[merge_list_used].ucMaxQS = usQS;
               (*lpSNVClustered)[merge_list_used].unCnt = 1;
               (*lpSNVClustered)[merge_list_used].bKeep = true;
               if (usPEnd != 0) (*lpSNVClustered)[merge_list_used].unPEndCnt = 1;
               if (cStrand==0) (*lpSNVClustered)[merge_list_used].unRvsCnt = 1;
               else if (cStrand==1) (*lpSNVClustered)[merge_list_used].unFwdCnt = 1;

               aulTotalRec[i]++;
          }
          else
          {
               (*lpSNVClustered)[merge_list_used].ulTotalQP += usQP;
               (*lpSNVClustered)[merge_list_used].ulTotalQS += usQS;

               if ((*lpSNVClustered)[merge_list_used].ucMaxQS < usQS)
                   (*lpSNVClustered)[merge_list_used].ucMaxQS = usQS;
            
               //if ((*lpSNVClustered)[merge_list_used].sMaxPEnd < usPEnd)
               //    (*lpSNVClustered)[merge_list_used].sMaxPEnd = usPEnd;

              (*lpSNVClustered)[merge_list_used].unCnt++;

               if (usPEnd != 0) (*lpSNVClustered)[merge_list_used].unPEndCnt++;

               if (cStrand==0) (*lpSNVClustered)[merge_list_used].unRvsCnt++;
               else if (cStrand==1) (*lpSNVClustered)[merge_list_used].unFwdCnt++;
          }
          if (++ulTest1 == 1000000){fprintf(stdout,"Total Recs gone through clustering = %luM\n", ++ulTest2); ulTest1=0;}
       }
    }

    fprintf(stdout,"Total Clustered Recs = %lu\n", merge_list_used+1);

    if (m_bOutClusteredFile) PrintClusteredRecs(lpSNVClustered,merge_list_used);
}


void CSXReadSNV::PrintClusteredRecs(stSNVClusterRec **lpSNVClustered, size_t merge_list_used)
{
    fprintf(stdout,"Print Clustered Recs ...\n");
    char acClusterFileName[1024]; sprintf(acClusterFileName,"%s/%s_clustered.txt",m_pcOutputDir,m_pcSampleID);
    FILE *pfCluster=fopen(acClusterFileName,"w");

    unsigned long ulTest1=0, ulTest2=0;  

    for (unsigned long i=0; i < merge_list_used+1; i++)
    {
         if (++ulTest1==100000) {fprintf(stdout,"Total Clustered Recs printed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}  

         fprintf(pfCluster,"%d\t%u\t%c>%c\t%u\t%lu\t%lu\t%u\n",
                 (*lpSNVClustered)[i].nChromosome,(*lpSNVClustered)[i].unOffset,
                 (*lpSNVClustered)[i].cRefBase,(*lpSNVClustered)[i].cVarBase,
                 (*lpSNVClustered)[i].unCnt,(*lpSNVClustered)[i].ulTotalQS,
                 (*lpSNVClustered)[i].ulTotalQP,(*lpSNVClustered)[i].ucMaxQS);
    }

    fclose(pfCluster);    
}


void CSXReadSNV::LogInvalidRec(FILE *pf,stSNVRec &Rec,int nChromosome)
{
    fprintf(pf,"%d\t%u\t%c\t%c\t%d\t%d\n",nChromosome, Rec.unOffset,
            Rec.cRefBase, Rec.cVarBase,Rec.ucQryPos,Rec.ucQltyScore);
}


void CSXReadSNV::FilterSNVTbl(size_t merge_list_used)
{
    for (unsigned long i=0; i < merge_list_used; i++)
    {
        if (m_pSNVClustered[i].ucMaxQS < m_nMinQScore )  //Ilumina - 30, CG-0
            m_pSNVClustered[i].bKeep=false;

        //if (m_pSNVClustered[i].sMaxPEnd < 1)
        //    m_pSNVClustered[i].bKeep=false;
    }

    unsigned long j,k;

    for (unsigned long i=0; i < merge_list_used; i++)
    {
        if (!m_pSNVClustered[i].bKeep) continue;

        j=i+1;

        if (!m_pSNVClustered[j].bKeep) continue;

        if (m_pSNVClustered[j].unOffset - m_pSNVClustered[i].unOffset > 1)
            continue;

        if (m_pSNVClustered[i].nChromosome!=m_pSNVClustered[j].nChromosome)
            continue;

        k=i+2;

        if (m_pSNVClustered[k].unOffset - m_pSNVClustered[j].unOffset > 1)
        {
            FilterSNVRec(m_pSNVClustered[i], m_pSNVClustered[j]); //comparing 2 set of data

        }
        else
        {
            FilterSNVRec(m_pSNVClustered[i], m_pSNVClustered[j],
                         m_pSNVClustered[k]);                     //comaparing 3 set of data

        }
    }
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

       fprintf(stdout,"Human Ref. File %d => ",i+1);
       for (int j=0; j < 70; j++)
       {
            fprintf(stdout,"%c",m_apcHsRef[i][j]);
       }
       fprintf(stdout,"...\n");
   }   
   CloseHsRefFiles();
}


bool CSXReadSNV::SearchClusteredList(unsigned long ulStart,unsigned long ulEnd,
                                     unsigned int unOffset, unsigned long &ulFound)
{
    unsigned long ulSearch = 0;

    while(!(ulEnd < ulStart)){
       ulSearch = (ulStart+ulEnd)/2;
       if (m_pSNVClustered[ulSearch].unOffset == unOffset){           
           ulFound=ulSearch; return true;
       }
       else if (m_pSNVClustered[ulSearch].unOffset < unOffset)
            ulStart = ulSearch+1;
       else{
            if (ulSearch > 0) ulEnd = ulSearch-1;
            else break;
       }

    };

    return false;
}


void CSXReadSNV::ProcAmbiguityList()
{
    char acFileName[1024]; sprintf(acFileName,"Ambiguity_%s.lst",m_pcSampleID);
    FILE *pfList = fopen(acFileName,"w");

    sprintf(acFileName,"Ambiguity_merge_%s.lst",m_pcSampleID);
    FILE *pfMergedList = fopen(acFileName,"w");

    int fd=open(m_pcAmbiguityFile, O_RDONLY, S_IRUSR);
    if (fd == -1) {fprintf(stderr, "File %s not found...\n",m_pcAmbiguityFile); return;}

    struct stat st; fstat(fd, &st);

    size_t file_size = (size_t)st.st_size; unsigned long ulTotalRecs = file_size/sizeof(stSNV);

    stSNV* pAmbiguitySNV=(stSNV*)mmap(NULL,file_size,PROT_READ,MAP_PRIVATE,fd,0);

    if (pAmbiguitySNV == (void *)-1) {
        fprintf(stderr, "Error: mapped %s failed: %s\n", m_pcAmbiguityFile,strerror(errno));return;
    }

    stSNVClusterRec *pAmbiguityClustered=NULL;unsigned long aulTotalRec[24];
    size_t merge_list_alloc=0, merge_list_used=0;

    GenerateSNVTblNew(pAmbiguitySNV,ulTotalRecs,file_size,&pAmbiguityClustered,aulTotalRec,
                      merge_list_alloc,merge_list_used);

    fprintf(stdout,"Merging Ambiguity Data...\n");

    unsigned long ulFound,ulTest1=0,ulTest2=0; bool bMerged = false;
    unsigned long ulIdx=0,j=0,ulCnt,ulStart=0,ulEnd,ulFoundSNP=0,ulFoundINS=0,ulFoundDEL=0;

    for (int i=0; i<24; i++)
    {
        if (i==0){ulEnd=m_aulTotalRec[i]-1;}
        else{ ulStart+=m_aulTotalRec[i-1]; ulEnd = ulStart+m_aulTotalRec[i]-1; }

        ulCnt =0;
        for (; ulCnt<aulTotalRec[i]; ulCnt++)
        {
            if (++ulTest1 > 1000000){
               ulTest1=0; fprintf(stdout,"Total Data checked for merging = %luM\n",++ulTest2);
            }

            bMerged = false;
            if (SearchClusteredList(ulStart,ulEnd,pAmbiguityClustered[ulIdx].unOffset,ulFound))
            {
                  for (j=ulFound-10; j < ulFound+10; j++)
                  {
                       if (m_pSNVClustered[j].unOffset > pAmbiguityClustered[ulIdx].unOffset) break;

                       if (m_pSNVClustered[j].unOffset == pAmbiguityClustered[ulIdx].unOffset &&
                           m_pSNVClustered[j].cRefBase == pAmbiguityClustered[ulIdx].cRefBase &&
                           m_pSNVClustered[j].cVarBase == pAmbiguityClustered[ulIdx].cVarBase
                       )
                       {
                           m_pSNVClustered[j].ulTotalQP += pAmbiguityClustered[ulIdx].ulTotalQP;
                           m_pSNVClustered[j].ulTotalQS += pAmbiguityClustered[ulIdx].ulTotalQS;

                           if (m_pSNVClustered[j].ucMaxQS < pAmbiguityClustered[ulIdx].ucMaxQS)
                               m_pSNVClustered[j].ucMaxQS = pAmbiguityClustered[ulIdx].ucMaxQS;

                           m_pSNVClustered[j].unCnt += pAmbiguityClustered[ulIdx].unCnt;

                           m_pSNVClustered[j].unPEndCnt += pAmbiguityClustered[ulIdx].unPEndCnt;

                           m_pSNVClustered[j].unRvsCnt += pAmbiguityClustered[ulIdx].unRvsCnt;
                           m_pSNVClustered[j].unFwdCnt += pAmbiguityClustered[ulIdx].unFwdCnt;

                           bMerged = true;
                           
                           if (m_pSNVClustered[j].cRefBase!='-' && m_pSNVClustered[j].cVarBase!='-') ulFoundSNP++;
                           else if (m_pSNVClustered[j].cRefBase!='-') ulFoundINS++;                           
                           else ulFoundDEL++;
                                                         
                           fprintf(pfMergedList,"%d\t%u\t%c>%c\t%u\t%u\n",
                           pAmbiguityClustered[ulIdx].nChromosome,pAmbiguityClustered[ulIdx].unOffset,
                           pAmbiguityClustered[ulIdx].cRefBase,pAmbiguityClustered[ulIdx].cVarBase,
                           m_pSNVClustered[j].unCnt,pAmbiguityClustered[ulIdx].unCnt);
                           
                           break;
                       }
                 }
            }
            if (!bMerged){
                fprintf(pfList,"%d\t%u\t%c>%c\t%u\t%lu\t%lu\t%u\n",
                pAmbiguityClustered[ulIdx].nChromosome,pAmbiguityClustered[ulIdx].unOffset,
                pAmbiguityClustered[ulIdx].cRefBase,pAmbiguityClustered[ulIdx].cVarBase,
                pAmbiguityClustered[ulIdx].unCnt,pAmbiguityClustered[ulIdx].ulTotalQS,
                pAmbiguityClustered[ulIdx].ulTotalQP,pAmbiguityClustered[ulIdx].ucMaxQS);
            }
            ulIdx++;
        }
    }

    fprintf(pfMergedList,
           "Merged SNP Counts = %lu\nMerged INS Counts = %lu\nMerged DEL Counts = %lu\n",            
           ulFoundSNP,ulFoundINS,ulFoundDEL);

    if (pAmbiguityClustered) free(pAmbiguityClustered);

    fclose(pfList); fclose(pfMergedList);

}


void CSXReadSNV::OutputAvgQScoreTable() 
{
    FILE *pfdbSnp=NULL, *pfdbIndel=NULL, *pfdbINS=NULL, *pfdbDEL=NULL; 

    char acFile[strlen(m_pcDensFilePath)+16], acFileName[1024];
    m_apfRD = new FILE*[24]; m_apfPDens = new FILE*[24];
    for (int i=0; i<24; i++){ m_apfRD[i]=NULL; m_apfPDens[i]=NULL; }     

    if (m_bWAmbiguityFiles){
        FILE * pf = fopen(m_pcAmbiguityFile,"r"); //Just to make sure ambiguity file exist.
        if (!pf){ fprintf(stdout,"Failed to open %s ...\n",m_pcAmbiguityFile); goto Exit;}
        fclose(pf); 
    }  
 
    pfdbSnp = fopen(m_pcdbSnp,"r"); if (!pfdbSnp){fprintf(stdout,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;}

    if (m_pcdbIndel){ 
        pfdbIndel = fopen(m_pcdbIndel,"r"); if (!pfdbSnp){fprintf(stdout,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;}
    }
    else {
       pfdbINS = fopen(m_pcdbINS,"r"); if (!pfdbINS) {fprintf(stderr,"Failed to open %s ...\n",m_pcdbINS); goto Exit;}
       pfdbDEL = fopen(m_pcdbDEL,"r"); if (!pfdbDEL) {fprintf(stderr,"Failed to open %s ...\n",m_pcdbDEL); goto Exit;}
    } 
 
    for (int i=0; i<24;i++)
    {
        sprintf(acFile,"%s/%s.seq%d.read_den_all",m_pcDensFilePath,m_pcProjSampleID,i+1);
        m_apfRD[i] = fopen(acFile,"r");
        if (!m_apfRD[i]){fprintf(stdout,"Failed to open %s ...\n", acFile); goto Exit;}
        else fprintf(stdout,"Touch %s ...\n",acFile);
    }

    for (int i=0; i<24;i++)
    {
        sprintf(acFile,"%s/%s.seq%d.read_den_perfect",m_pcDensFilePath,m_pcProjSampleID,i+1);
        m_apfPDens[i] = fopen(acFile,"r");
        if (!m_apfPDens[i]){fprintf(stdout,"Failed to open %s ...\n", acFile); goto Exit;}
        else fprintf(stdout,"Touch %s ...\n",acFile);
    }

    m_apfHsRef = new FILE*[24]; m_apcHsRef = new char*[24];

    for (short i=0;i<24;i++)
    {
        if (i < 22){sprintf(acFileName,"%s/hs_ref_GRCh37_chr%d.fa",m_pcHsRefFilePath,i+1);}
        else if(i == 22) sprintf(acFileName,"%s/hs_ref_GRCh37_chrX.fa",m_pcHsRefFilePath);
        else sprintf(acFileName,"%s/hs_ref_GRCh37_chrY.fa",m_pcHsRefFilePath);

        m_apfHsRef[i] = fopen(acFileName,"r");

        if (!m_apfHsRef[i]) {fprintf(stdout,"Failed to open %s ....\n",acFileName); goto Exit;}
        else fprintf(stdout,"Touch %s ...\n",acFileName);        
    }

    fprintf(stdout,"Generating %s Map...\n",m_pcdbSnp);
    
    if (m_pcdbIndel) {
        if (!GenerateDBSNPMap_OdbSnpfmt(pfdbSnp)) goto Exit;

        fprintf(stdout,"Generating %s Map...\n",m_pcdbIndel);
        if (!GenerateDBINDELMap(pfdbIndel)) goto Exit;
    }
    else{
        if (!GenerateDBSNPMap(pfdbSnp)) goto Exit;

        fprintf(stdout,"Generating %s Map...\n",m_pcdbINS);
        m_pObjdbINS = new CSXDBIndelgChecker(); if (!m_pObjdbINS) goto Exit;          
        if (!m_pObjdbINS->GenerateDBINDELMap(pfdbINS)) goto Exit;
           
        fprintf(stdout,"Generating %s Map...\n",m_pcdbDEL);
        m_pObjdbDEL = new CSXDBIndelgChecker(); if (!m_pObjdbDEL) goto Exit;
        if (!m_pObjdbDEL->GenerateDBINDELMap(pfdbDEL)) goto Exit;
    }

    fprintf(stdout,"Generating SNV Table...\n"); 

    for (int i=0; i<24;i++) {m_aulTotalRec[i]=0;}

    GenerateSNVTblNew(m_pSNV,m_ulTotalRecs,m_file_size,&m_pSNVClustered,m_aulTotalRec,
                      m_merge_list_alloc,m_merge_list_used);

    if (m_bWAmbiguityFiles) {fprintf(stdout,"Process Ambiguity List....\n"); ProcAmbiguityList();}

    //fprintf(stdout,"Filtering SNV Table...\n"); FilterSNVTbl(m_merge_list_used);
    fprintf(stdout,"Load Human Reference Sequence ...\n"); LoadHRefSeq();
    fprintf(stdout,"Printing SNP, INS, DEL Lists...\n"); PrintAvgQScoreList();

Exit:

    if (pfdbSnp) fclose(pfdbSnp); if (pfdbIndel) fclose(pfdbIndel);
    ClrDBSNPMap(); ClrDBINDELMap();

    for (int j=0;j<24; j++){
         if (m_apfRD[j]) fclose(m_apfRD[j]);
         if (m_apfPDens[j]) fclose(m_apfPDens[j]);
    }

    delete[] m_apfRD; delete[] m_apfPDens;
}


int CSXReadSNV::CheckPoly(int nChroIdx, unsigned int unOffset, char cAllele,eSNVType eSNV,FILE *pf)
{
    //fprintf(pf,"%d\t%d\t=>\t",nChroIdx+1,unOffset);

    if (eSNV == eDEL) unOffset--; // if it is an INS we need to shift 1 bp to the right. Thus, no need to minus the offset by 1
       
    if (m_apcHsRef[nChroIdx][unOffset]!=cAllele) { 
       if (eSNV == eDEL)unOffset++; 
      
       //fprintf(pf,"%d\t%c\n",unOffset,cAllele);
       return unOffset;
    }

    unsigned int i;   

    for (i=unOffset+1;i < m_anHsRefSize[nChroIdx]; i++)
    {
        if (m_apcHsRef[nChroIdx][i] != cAllele)
           break; 
    }

    /////Testing///////////////////////////////
    /*
    if (unOffset != i-1){
       int nSize = i - (unOffset+1); 
       char acAllele[nSize+1];
       memcpy(&acAllele[0],&m_apcHsRef[nChroIdx][unOffset],nSize);
       acAllele[nSize]=0;
       
       //fprintf(pf,"%d\t%c\t%s\n",i,cAllele,acAllele);      
    }
    else
    {
       fprintf(pf,"*%d\t%c\t%c\n",i,cAllele,m_apcHsRef[nChroIdx][unOffset]);
    }*/  

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
    char acPDensFName[1024]; sprintf(acPDensFName,"%s/PDens_del_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenDEL = fopen(acPDensFName,"w");

    sprintf(acPDensFName,"%s/PDens_ins_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenINS = fopen(acPDensFName,"w");

    stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astDELStats[m_unMaxSupportingReads+1];
    int nSNPIdx,nINSIdx,nDELIdx=0;

    unsigned long ulTest1=0, ulTest2=0;    
    
    unsigned int unChroIdx=0,unReadDensity=0,unTotPDens=0;    
    bool bdbSNP=false; char acRsid[4096];    
    float fAvgQryPos,fAvgQScore;
    
    FILE *pfSNP,*pfINS,*pfDEL,*pfSyno,*pfSNPStats,*pfINSStats,*pfDELStats;
    pfSNP=pfINS=pfDEL=pfSyno=pfSNPStats=0,pfINSStats=pfDELStats=0;

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[1024];

    sprintf(acFileName,"%s/%s_snp_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}


    sprintf(acFileName,"%s/%s_snp_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;} 
   
///////

    sprintf(acFileName,"%s/%s_ins_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_ins_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////
    
    sprintf(acFileName,"%s/%s_del_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_del_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;} 

///////

    fprintf(pfSNP,"#REPORT NAME\tSNP Avg. Score Distribution\n");
    fprintf(pfSNP,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfSNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfSNP,"#LANE NO\tALL\n");
    fprintf(pfSNP,"#GENERATED AT\t"); PrintRptDateTime(pfSNP);
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfSNP,"#FILTER\t>=1\n\n\n");

    fprintf(pfSNP,"Chromosome\tOffset\tNucleotide_Variant\tFwd_SNP_Reads"
                  "\tRvs_SNP_Reads\tTotal_SNP_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");

    fprintf(pfINS,"#REPORT NAME\tINS Avg. Score Distribution\n");
    fprintf(pfINS,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
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
    fprintf(pfDEL,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfDEL,"Chromosome\tOffset\tDeleted_Base\tFwd_DEL_Reads"
                  "\tRvs_DEL_Reads\tTotal_DEL_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n"); 

    for (unsigned long i=0; i < (m_merge_list_used+1); i++)
    {
        if (++ulTest1 == 1000000){fprintf(stdout,"Total Clustered Recs Processed = %luM\n", ++ulTest2); ulTest1=0; }

        if (m_pSNVClustered[i].nChromosome == 0) continue;

        unChroIdx = m_pSNVClustered[i].nChromosome-1;

        if (m_pSNVClustered[i].cRefBase!='-' && m_pSNVClustered[i].cVarBase!='-')
        {
            bdbSNP=OutputSNPNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                  m_pSNVClustered[i].cVarBase,pfSNP,acRsid);

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfSNP,m_pSNVClustered[i].nChromosome);
            fprintf(pfSNP,"\t%u\t%c>%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,
                    m_pSNVClustered[i].unCnt);

            fprintf(pfSNP,"\t%u\t%u\t%.2f\t%.2f\t%s",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfSNP,"\n");

            nSNPIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astSNPStats[nSNPIdx].unCnt++;
            if (bdbSNP) astSNPStats[nSNPIdx].undbSNP++;
            astSNPStats[nSNPIdx].ulTotalReadDens += unReadDensity;
        }
        else if (m_pSNVClustered[i].cRefBase =='-' )
        {
            if (m_pcdbIndel) bdbSNP=OutputINDELNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                                     m_pSNVClustered[i].cVarBase, pfINS,acRsid);                                                    else bdbSNP=m_pObjdbINS->OutputINDELNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                                      m_pSNVClustered[i].cVarBase, pfINS,acRsid); 

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;
            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;

            unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx],m_pSNVClustered[i].unOffset);

            if ((Percentage(m_pSNVClustered[i].unCnt,unReadDensity)<m_fMinReadStrength)&& !bdbSNP) continue;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);

            OutputChrom(pfINS,m_pSNVClustered[i].nChromosome);

            fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,
                    m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,m_pSNVClustered[i].unCnt);

            fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfINS,"\n");

            nINSIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;

            astINSStats[nINSIdx].unCnt++; if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
            astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;

            OutputChrom(pfPDenINS,m_pSNVClustered[i].nChromosome);

            unTotPDens = GetReadDensityEx(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);

            fprintf(pfPDenINS,"\t%u\t%u\t-\t%u\t%u\n",m_pSNVClustered[i].unOffset,
                              unTotPDens+m_pSNVClustered[i].unCnt,unTotPDens,
                              m_pSNVClustered[i].unCnt);
        } 
        else
        {
            if (m_pcdbIndel) bdbSNP=OutputINDELNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                                     m_pSNVClustered[i].cVarBase, pfDEL,acRsid);
            else bdbSNP=m_pObjdbDEL->OutputINDELNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                                      m_pSNVClustered[i].cVarBase, pfDEL,acRsid);

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;
            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;

            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset) + m_pSNVClustered[i].unCnt;        

            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;
            
            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);

            OutputChrom(pfDEL,m_pSNVClustered[i].nChromosome);
            fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,m_pSNVClustered[i].unCnt);

            fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s",
                    unReadDensity,m_pSNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfDEL,"\n");

            nDELIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;

            astDELStats[nDELIdx].unCnt++; if (bdbSNP) astDELStats[nDELIdx].undbSNP++;
            astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;

            OutputChrom(pfPDenDEL,m_pSNVClustered[i].nChromosome);
            unTotPDens = GetReadDensity(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);

            fprintf(pfPDenDEL,"\t%u\t%u\t-\t%u\t%u\n",m_pSNVClustered[i].unOffset,
                              unTotPDens+m_pSNVClustered[i].unCnt,
                              unTotPDens,m_pSNVClustered[i].unCnt); 
        }
    }//end for

//////////////////////

    fprintf(stderr,"Printing AvgQScore Statistic Report...\n");

    PrintAvgQSStatsRpt(pfSNPStats,astSNPStats,(char*)"SNP");
    PrintAvgQSStatsRpt(pfINSStats,astINSStats,(char*)"INS");
    PrintAvgQSStatsRpt(pfDELStats,astDELStats,(char*)"DEL");

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfINS) fclose(pfINS);
    if (pfDEL) fclose(pfDEL); if (pfSNPStats) fclose(pfSNPStats);
    if (pfINSStats) fclose(pfINSStats); if (pfDELStats) fclose(pfDELStats);
    if (pfPDenINS) fclose(pfPDenINS); if (pfPDenDEL) fclose(pfPDenDEL); 
}

/*
void CSXReadSNV::PrintAvgQScoreList()
{
    char acPolyINSFileName[1024]; sprintf(acPolyINSFileName,"%s/PolyINS_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPolyINS = fopen(acPolyINSFileName,"w");

    char acPolyDELFileName[1024]; sprintf(acPolyDELFileName,"%s/PolyDEL_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPolyDEL = fopen(acPolyDELFileName,"w"); 

    char acPDensFName[1024]; sprintf(acPDensFName,"%s/PDens_del_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenDEL = fopen(acPDensFName,"w");

    sprintf(acPDensFName,"%s/PDens_ins_%s.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenINS = fopen(acPDensFName,"w");

    stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astDELStats[m_unMaxSupportingReads+1];
    int nSNPIdx,nINSIdx,nDELIdx=0;   

    unsigned long ulTest1=0, ulTest2=0;
    stChkPoly ChkPoly_INS, ChkPoly_DEL;
    bool bINS1=true, bDEL1=true;
    unsigned int unChroIdx=0,unReadDensity=0,unOffset=0,unTotPDens=0;

    bool bdbSNP=false; char acRsid[4096];

    float fAvgQryPos,fAvgQScore;

    FILE *pfSNP,*pfINS,*pfDEL,*pfSyno,*pfSNPStats,*pfINSStats,*pfDELStats;
    pfSNP=pfINS=pfDEL=pfSyno=pfSNPStats=0,pfINSStats=pfDELStats=0;

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[24+strlen(m_pcSampleID)];
    sprintf(acFileName,"%s/%s_snp_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNP = fopen(acFileName,"w");
    if (!pfSNP) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}


    sprintf(acFileName,"%s/%s_snp_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfSNPStats = fopen(acFileName,"w");
    if (!pfSNPStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    sprintf(acFileName,"%s/%s_ins_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINS = fopen(acFileName,"w");
    if (!pfINS) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_ins_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfINSStats = fopen(acFileName,"w");
    if (!pfINSStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////
    sprintf(acFileName,"%s/%s_del_c1_aqs_%02d%02d%02d.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDEL = fopen(acFileName,"w");
    if (!pfDEL) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_del_summary_c1_aqs_%02d%02d%02d.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfDELStats = fopen(acFileName,"w");
    if (!pfDELStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    fprintf(pfSNP,"#REPORT NAME\tSNP Avg. Score Distribution\n");
    fprintf(pfSNP,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
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
    fprintf(pfINS,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
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
    fprintf(pfDEL,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfDEL,"Chromosome\tOffset\tDeleted_Base\tFwd_DEL_Reads"                  
                  "\tRvs_DEL_Reads\tTotal_DEL_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");
 
    for (unsigned long i=0; i < (m_merge_list_used+1); i++)
    {        
        if (++ulTest1 == 1000000){fprintf(stdout,"Total Clustered Recs Processed = %luM\n", ++ulTest2); ulTest1=0; }

        if (m_pSNVClustered[i].nChromosome == 0) continue;

        unChroIdx = m_pSNVClustered[i].nChromosome-1;    

        if (m_pSNVClustered[i].cRefBase!='-' && m_pSNVClustered[i].cVarBase!='-')
        {
            bdbSNP=OutputSNPNovel(m_pSNVClustered[i].nChromosome,m_pSNVClustered[i].unOffset,
                                  m_pSNVClustered[i].cVarBase,pfSNP,acRsid);

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[i].ulTotalQP,m_pSNVClustered[i].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);            
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue; 
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfSNP,m_pSNVClustered[i].nChromosome);
            fprintf(pfSNP,"\t%u\t%c>%c\t%u\t%u\t%u",m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,
                    m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unFwdCnt,m_pSNVClustered[i].unRvsCnt,
                    m_pSNVClustered[i].unCnt);

            fprintf(pfSNP,"\t%u\t%u\t%.2f\t%.2f\t%s",
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
                bINS1=false;ChkPoly_INS.ulIdx=i; ChkPoly_INS.bPoly=false;
                ChkPoly_INS.unOffset=CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,eINS,pfPolyINS);

                if (ChkPoly_INS.unOffset!=m_pSNVClustered[i].unOffset){ m_pSNVClustered[i].bKeep = true; ChkPoly_INS.bPoly=true;}

                ChkPoly_INS.unPDensity = GetReadDensityEx(m_apfPDens[unChroIdx],m_pSNVClustered[i].unOffset);
                ChkPoly_INS.unCnt = 1;
                if (ChkPoly_INS.unPDensity > 5)ChkPoly_INS.unPDensity -=5;  

                continue;
            }

            unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cVarBase,eINS,pfPolyINS);            

            if (m_pSNVClustered[i].nChromosome==m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome && 
                unOffset == ChkPoly_INS.unOffset && 
                m_pSNVClustered[i].cVarBase==m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase)
            { 
               if (ChkPoly_INS.ulIdx !=i)
               { 
                  AddReadCnt(i,ChkPoly_INS.ulIdx); ChkPoly_INS.ulIdx=i;

                  unReadDensity = GetReadDensityEx(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);
                  ChkPoly_INS.unCnt++;
                  if (unReadDensity > 5) unReadDensity -=5;
                  ChkPoly_INS.unPDensity += unReadDensity; 
               }

               m_pSNVClustered[ChkPoly_INS.ulIdx].bKeep = true; ChkPoly_INS.bPoly=true;
               continue; 
            }
            
            if (m_pcdbIndel) bdbSNP=OutputINDELNovel(m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome,ChkPoly_INS.unOffset,
                                                      m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase, pfINS,acRsid);                                          else bdbSNP=m_pObjdbINS->OutputINDELNovel(m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome,ChkPoly_INS.unOffset,
                                                      m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase, pfINS,acRsid);

            if (!m_pSNVClustered[ChkPoly_INS.ulIdx].bKeep && !bdbSNP) goto EndINS;

            if (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP ) {
                goto EndINS;
            } 

            unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx],ChkPoly_INS.unOffset);

            if ((Percentage(m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt,unReadDensity)<m_fMinReadStrength)&& !bdbSNP) 
               goto EndINS;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_INS.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

            OutputChrom(pfINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);

            fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",ChkPoly_INS.unOffset,m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase,
                    m_pSNVClustered[ChkPoly_INS.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_INS.ulIdx].unRvsCnt,
                    m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

            fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s",
                    unReadDensity,m_pSNVClustered[ChkPoly_INS.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);
            
            fprintf(pfINS,"\n");

            nINSIdx = (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt-1;

            astINSStats[nINSIdx].unCnt++; if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
            astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;
            
            if (ChkPoly_INS.bPoly){
               OutputChrom(pfPolyINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);
               fprintf(pfPolyINS,"\t%u\n",ChkPoly_INS.unOffset);
            }
            
            OutputChrom(pfPDenINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);
 
            if (ChkPoly_INS.unPDensity > 0) unTotPDens = int(((float)ChkPoly_INS.unPDensity/(float)ChkPoly_INS.unCnt)+0.5);
            else unTotPDens = ChkPoly_INS.unPDensity;

            fprintf(pfPDenINS,"\t%u\t%u\t-\t%u\t%u\t%u\n",ChkPoly_INS.unOffset,
                              unTotPDens+m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt,
                              ChkPoly_INS.unPDensity,
                              ChkPoly_INS.unCnt,
                              m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);
EndINS:
            ChkPoly_INS.ulIdx = i; ChkPoly_INS.unOffset = unOffset; 

            ChkPoly_INS.unPDensity = GetReadDensityEx(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);
            ChkPoly_INS.unCnt = 1;

            if (ChkPoly_INS.unPDensity > 5) ChkPoly_INS.unPDensity -=5; 
  
            if (ChkPoly_INS.unOffset!=m_pSNVClustered[i].unOffset) {m_pSNVClustered[i].bKeep=true; ChkPoly_INS.bPoly=true;}
            else ChkPoly_INS.bPoly=false;    
        }
        else
        {
	    if (bDEL1){
                bDEL1=false;ChkPoly_DEL.ulIdx=i; ChkPoly_DEL.bPoly=false;
                ChkPoly_DEL.unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,eDEL,pfPolyDEL);
                ChkPoly_DEL.unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset); 

                if (ChkPoly_DEL.unOffset != m_pSNVClustered[i].unOffset) 
                {
                    m_pSNVClustered[i].bKeep = true; ChkPoly_DEL.bPoly=true;
 
                    if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt) ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;
                    else ChkPoly_DEL.unReadDensity = 0;
                } 

                ChkPoly_DEL.unPDensity = GetReadDensity(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);
                ChkPoly_DEL.unCnt = 1;
                if (ChkPoly_DEL.unPDensity > 5)ChkPoly_DEL.unPDensity -=5;
 
                continue;
            }

            unOffset = CheckPoly(unChroIdx,m_pSNVClustered[i].unOffset,m_pSNVClustered[i].cRefBase,eDEL,pfPolyDEL);    
            
            if (m_pSNVClustered[i].nChromosome==m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome && 
                unOffset == ChkPoly_DEL.unOffset &&
                m_pSNVClustered[i].cRefBase==m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase)
            {
               if (ChkPoly_DEL.ulIdx!=i) 
               {
                   AddReadCnt(i,ChkPoly_DEL.ulIdx);ChkPoly_DEL.ulIdx=i;

                   unReadDensity = GetReadDensity(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);
                   ChkPoly_DEL.unCnt++;
                   if (unReadDensity > 5) unReadDensity -=5;
                   ChkPoly_DEL.unPDensity += unReadDensity; 
               }           
               
               if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt) ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;
               else ChkPoly_DEL.unReadDensity = 0;
               
               m_pSNVClustered[i].bKeep = true; ChkPoly_DEL.bPoly=true;
               continue;
            }

            if (m_pcdbIndel) bdbSNP=OutputINDELNovel(m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome,ChkPoly_DEL.unOffset,
                                                     m_pSNVClustered[ChkPoly_DEL.ulIdx].cVarBase, pfDEL,acRsid); 
            else bdbSNP=m_pObjdbDEL->OutputINDELNovel(m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome,ChkPoly_DEL.unOffset,
                                                      m_pSNVClustered[ChkPoly_DEL.ulIdx].cVarBase, pfDEL,acRsid);
 
            if (!m_pSNVClustered[ChkPoly_DEL.ulIdx].bKeep && !bdbSNP) goto EndDEL;

            if (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP) goto EndDEL;

	    unReadDensity = ChkPoly_DEL.unReadDensity + m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt;

            if (Percentage(m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) goto EndDEL;

            fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
            fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

            OutputChrom(pfDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);
            fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",ChkPoly_DEL.unOffset,m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase,
                    m_pSNVClustered[ChkPoly_DEL.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_DEL.ulIdx].unRvsCnt,
                    m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
            
            fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s",
                    unReadDensity,m_pSNVClustered[ChkPoly_DEL.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

            fprintf(pfDEL,"\n");

            nDELIdx = (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads)?
                       m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt-1;

            astDELStats[nDELIdx].unCnt++; if (bdbSNP) astDELStats[nDELIdx].undbSNP++; 
            astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;

            if (ChkPoly_DEL.bPoly){
               OutputChrom(pfPolyDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);
               fprintf(pfPolyDEL,"\t%u\n",ChkPoly_DEL.unOffset);
            }

            OutputChrom(pfPDenDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);

            if (ChkPoly_DEL.unPDensity > 0) unTotPDens = int(((float)ChkPoly_DEL.unPDensity/(float)ChkPoly_DEL.unCnt)+0.5);
            else unTotPDens = ChkPoly_DEL.unPDensity;
    
            fprintf(pfPDenDEL,"\t%u\t%u\t-\t%u\t%u\t%u\n",ChkPoly_DEL.unOffset,
                              unTotPDens+m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,
                              ChkPoly_DEL.unPDensity,
                              ChkPoly_DEL.unCnt,
                              m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
EndDEL:
            ChkPoly_DEL.ulIdx = i; ChkPoly_DEL.unOffset = unOffset; 

            if (ChkPoly_DEL.unOffset!=m_pSNVClustered[i].unOffset) {m_pSNVClustered[i].bKeep=true; ChkPoly_DEL.bPoly=true;}
            else ChkPoly_DEL.bPoly=false;

            ChkPoly_DEL.unReadDensity = GetReadDensity(m_apfRD[unChroIdx],unOffset);
             
            if (ChkPoly_DEL.unOffset != m_pSNVClustered[i].unOffset)
            {
               if (ChkPoly_DEL.unReadDensity > m_pSNVClustered[i].unCnt){ 
                   ChkPoly_DEL.unReadDensity -= m_pSNVClustered[i].unCnt;  m_pSNVClustered[i].bKeep = true;}
               else ChkPoly_DEL.unReadDensity = 0;
            }             

            ChkPoly_DEL.unPDensity = GetReadDensity(m_apfPDens[unChroIdx], m_pSNVClustered[i].unOffset);
            ChkPoly_DEL.unCnt = 1;
            if (ChkPoly_DEL.unPDensity > 5)ChkPoly_DEL.unPDensity -=5;
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

    OutputChrom(pfINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);

    fprintf(pfINS,"\t%u\t%c\t%u\t%u\t%u",
                   m_pSNVClustered[ChkPoly_INS.ulIdx].unOffset,m_pSNVClustered[ChkPoly_INS.ulIdx].cVarBase,
                   m_pSNVClustered[ChkPoly_INS.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_INS.ulIdx].unRvsCnt,
                   m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);

    fprintf(pfINS,"\t%u\t%u\t%.2f\t%.2f\t%s",
                   unReadDensity,m_pSNVClustered[ChkPoly_INS.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

    fprintf(pfINS,"\n");

    nINSIdx = (m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt > m_unMaxSupportingReads)?
               m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt-1;

    astINSStats[nINSIdx].unCnt++;
    if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
    astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;

    if (ChkPoly_INS.bPoly){
               OutputChrom(pfPolyINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);
               fprintf(pfPolyINS,"\t%u\n",ChkPoly_INS.unOffset);
    }

    OutputChrom(pfPDenINS,m_pSNVClustered[ChkPoly_INS.ulIdx].nChromosome);

    if (ChkPoly_INS.unPDensity > 0) unTotPDens = int(((float)ChkPoly_INS.unPDensity/(float)ChkPoly_INS.unCnt)+0.5);
    else unTotPDens = ChkPoly_INS.unPDensity;

    fprintf(pfPDenINS,"\t%u\t%u\t-\t%u\t%u\t%u\n",ChkPoly_INS.unOffset,
                       unTotPDens+m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt,
                       ChkPoly_INS.unPDensity,ChkPoly_INS.unCnt,
                       m_pSNVClustered[ChkPoly_INS.ulIdx].unCnt);  

ProcLastDELRec:

    if (!m_pSNVClustered[ChkPoly_DEL.ulIdx].bKeep && !bdbSNP) goto PrintStatsRpt;

    fAvgQryPos = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQP,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);
    fAvgQScore = CalcAvgQS(m_pSNVClustered[ChkPoly_DEL.ulIdx].ulTotalQS,m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

    unReadDensity = ChkPoly_DEL.unReadDensity + m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt; 

    if (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads && !bdbSNP) goto PrintStatsRpt;
    if (Percentage(m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) goto PrintStatsRpt;

    OutputChrom(pfDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);
    fprintf(pfDEL,"\t%u\t%c\t%u\t%u\t%u",
                   m_pSNVClustered[ChkPoly_DEL.ulIdx].unOffset,m_pSNVClustered[ChkPoly_DEL.ulIdx].cRefBase,
                   m_pSNVClustered[ChkPoly_DEL.ulIdx].unFwdCnt,m_pSNVClustered[ChkPoly_DEL.ulIdx].unRvsCnt,
                   m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);

    fprintf(pfDEL,"\t%u\t%u\t%.2f\t%.2f\t%s",
                   unReadDensity,m_pSNVClustered[ChkPoly_DEL.ulIdx].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

    fprintf(pfDEL,"\n");

    nDELIdx = (m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt > m_unMaxSupportingReads)?
               m_unMaxSupportingReads:m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt-1;

    astDELStats[nDELIdx].unCnt++;

    if (bdbSNP) astDELStats[nDELIdx].undbSNP++;
    astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;

    if (ChkPoly_DEL.bPoly){
        OutputChrom(pfPolyDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);
        fprintf(pfPolyDEL,"\t%u\n",ChkPoly_DEL.unOffset);
    }

    OutputChrom(pfPDenDEL,m_pSNVClustered[ChkPoly_DEL.ulIdx].nChromosome);

    if (ChkPoly_DEL.unPDensity > 0) unTotPDens = int(((float)ChkPoly_DEL.unPDensity/(float)ChkPoly_DEL.unCnt)+0.5);
    else unTotPDens = ChkPoly_DEL.unPDensity;

    fprintf(pfPDenDEL,"\t%u\t%u\t-\t%u\t%u\t%u\n",ChkPoly_DEL.unOffset,
                      unTotPDens+m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt,
                      ChkPoly_DEL.unPDensity,ChkPoly_DEL.unCnt,
                      m_pSNVClustered[ChkPoly_DEL.ulIdx].unCnt);  

//////////////////////

PrintStatsRpt:

   fprintf(stderr,"Printing AvgQScore Statistic Report...\n");

    PrintAvgQSStatsRpt(pfSNPStats,astSNPStats,(char*)"SNP");
    PrintAvgQSStatsRpt(pfINSStats,astINSStats,(char*)"INS");
    PrintAvgQSStatsRpt(pfDELStats,astDELStats,(char*)"DEL");

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfINS) fclose(pfINS);
    if (pfDEL) fclose(pfDEL); if (pfSNPStats) fclose(pfSNPStats);
    if (pfINSStats) fclose(pfINSStats); if (pfDELStats) fclose(pfDELStats);
    if (pfPolyINS) fclose(pfPolyINS); if (pfPolyDEL) fclose(pfPolyDEL);
    if (pfPDenINS) fclose(pfPDenINS); if (pfPDenDEL) fclose(pfPDenDEL);
}
*/


void CSXReadSNV::PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{   
    fprintf(pf,"#REPORT NAME\t%s Avg. QScore Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pf,"#DISPLAY RANGE\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
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
               ">=%u\t%lu\t%.2f\t%lu\t%.2f\n",
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
                   ">=%u\t%lu\t%.2f\t%lu\t%.2f\n",
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

    fprintf(pf,"Total\t%lu\t\t%lu\n",
            m_stSNVAnalysis.ulTotalSNP,m_stSNVAnalysis.ulTotalSNP_dbsnp);
    fprintf(pf,"\n\n\n");
}


extern "C" int CompareMNVRecs(const void *a, const void *b)
{
    const stMNVRec *pa = (const stMNVRec *)a, *pb = (const stMNVRec *)b;

    if (pa->unOffset != pb->unOffset) return (pa->unOffset - pb->unOffset);

    int n = strcmp(pa->pcRefBase,pb->pcRefBase);

    if (n != 0) return (n > 0);

    return (strcmp(pa->pcVarBase,pb->pcVarBase) > 0);
}


/*
extern "C" int CompareMNVRecs(const void *a, const void *b)
{
    const stMNVRec *pa = (const stMNVRec *)a, *pb = (const stMNVRec *)b;     
    
    if (pa->unOffset != pb->unOffset) return (pa->unOffset - pb->unOffset);


    int nBps_a=0,nBps_b=0,nStop_a=0,nStop_b=0; 
    
    if (pa->pcVarBase[0]=='-' && pb->pcVarBase[0]=='-'){   //DEL
        nBps_a = strlen(pa->pcRefBase); nBps_b = strlen(pb->pcRefBase);
    } 
    else {
        nBps_a = strlen(pa->pcVarBase); nBps_b = strlen(pb->pcVarBase);
    } 

    nStop_a = (pa->unOffset+nBps_a)-1;
    nStop_b = (pb->unOffset+nBps_b)-1; 

    if (nStop_a != nStop_b) return (nStop_a - nStop_b); 

    int n = strcmp(pa->pcRefBase,pb->pcRefBase); 

    if (n != 0) return (n < 0);

    return (strcmp(pa->pcVarBase,pb->pcVarBase) < 0);
}
*/

void *SortMNVData(void *ptr)
{
    char tmpChr; stMNVRecList *pMNVRecList = (stMNVRecList*)ptr; stMNVRec *pstMNVRec = NULL;
    

    for (int i=0; i < pMNVRecList->nUsed; i++)
    {        
        pstMNVRec = &(pMNVRecList->pMNVRec[i]);

        for (unsigned int j=0; j < strlen(pstMNVRec->pcRefBase); j++)  
        { 
             tmpChr = toupper(pstMNVRec->pcRefBase[j]);
             pstMNVRec->pcRefBase[j] = tmpChr;

             tmpChr = toupper(pstMNVRec->pcVarBase[j]);
             pstMNVRec->pcVarBase[j] = tmpChr;
       }  
    }

    qsort(pMNVRecList->pMNVRec,pMNVRecList->nUsed,sizeof(stMNVRec), CompareMNVRecs);
    return 0;
}


void CSXReadSNV::AllocMNVMergeListSize(stMNVClusterRec **pMNVClustered, size_t &merge_list_alloc, size_t &merge_list_used)
{
        if (*pMNVClustered == NULL) {
                merge_list_alloc = INIT_MERGED_LIST_SIZE;
                merge_list_used = 0;
                try {
                        *pMNVClustered = (stMNVClusterRec*)malloc(sizeof(stMNVClusterRec)*merge_list_alloc);
                }
                catch (...) {
                        fprintf(stderr, "Error: insufficient memory for clustered MNV records list!\n");
                        fprintf(stderr, "Program terminated.\n");
                        exit(1);
                }

                return;
        }

        merge_list_alloc += EXT_MERGED_LIST_SIZE;
        try {
                *pMNVClustered = (stMNVClusterRec*)realloc(*pMNVClustered,sizeof(stMNVClusterRec)*merge_list_alloc);
        }
        catch (...) {
                fprintf(stderr, "Error: insufficient memory for clustered MNV records list!\n");
                fprintf(stderr, "Program terminated.\n");
                exit(1);
        }
}


void CSXReadSNV::PrintUnSortedMNV(FILE *pfInFile)
{
    fprintf(stdout,"Printing Unsorted MNV recs...\n");

    char  *pcAllele=NULL, *pcRefBps=NULL, *pcVarBps=NULL, acFileName[1024], cStrand; 
    sprintf(acFileName,"%s/%s_unsorted_MNV.txt",m_pcOutputDir,m_pcSampleID);

    int nChro; 
    unsigned char ucCnt, ucQryPos, *pucQltyScore=NULL; 
    unsigned short usPEnd; unsigned int unOffset;
    FILE *pfIn=fopen(acFileName,"w");
    unsigned long ulFileOfs, ulTest1=0, ulTest2=0; 
    
    while (!(feof(pfInFile)))
    {
        nChro=0; fread(&nChro,4,1,pfInFile); if (nChro==0) break;
       
        if (++ulTest1==100000){fprintf(stdout,"Total Unsorted Recs being printed = %lu X 1.0M\n",++ulTest2);ulTest1=0;}

        fprintf(pfIn,"%d\t",nChro);                    

        fread(&unOffset,4,1,pfInFile); fprintf(pfIn,"%u\t",unOffset);              //fprintf(stdout,"Offset-%u\t",unOffset);
        fread(&ucCnt,1,1,pfInFile); fprintf(pfIn,"%d\t",ucCnt);                    //fprintf(stdout,"Count-%d\t",ucCnt);

        ulFileOfs = ftell(pfInFile);

        pcAllele=(char*)malloc(ucCnt+2);

        fread(pcAllele,ucCnt+1,1,pfInFile); pcAllele[ucCnt+1]=0;

        if (pcAllele[0]=='-')//{//INS
            fprintf(pfIn,"-\t%s\t",pcAllele+1);  //fprintf(stdout,"RefBase- -\tVarBase-%s\t",pcAllele+1);}
        else if (pcAllele[ucCnt]=='-') {//DEL
            pcAllele[ucCnt] = 0; fprintf(pfIn,"%s\t-\t",pcAllele);  //fprintf(stdout,"RefBase-%s\tVarBase- -\t",pcAllele);              
        } else { //SNP
            pcRefBps = (char*)malloc(ucCnt+1);
            pcVarBps = (char*)malloc(ucCnt+1); 

            fseek(pfInFile,ulFileOfs,SEEK_SET); 

            fread(pcRefBps,ucCnt,1,pfInFile); pcRefBps[ucCnt] = 0;
            fread(pcVarBps,ucCnt,1,pfInFile); pcVarBps[ucCnt] = 0;

            fprintf(pfIn,"%s\t%s\t",pcRefBps,pcVarBps);  //fprintf(stdout,"RefBase-%s\tVarBase-%s\t",pcRefBps,pcVarBps);
        }  
        
        if (pcAllele) {free(pcAllele); pcAllele=NULL;}
        if (pcRefBps) {free(pcRefBps); pcRefBps=NULL;}
        if (pcVarBps) {free(pcVarBps); pcVarBps=NULL;}

        fread(&ucQryPos,1,1,pfInFile); fprintf(pfIn,"%d\t",ucQryPos); //fprintf(stdout,"QryPos-%d\t",ucQryPos);
        
        pucQltyScore = (unsigned char*)malloc(1);

        //fread(pucQltyScore,ucCnt,1,pfInFile); pucQltyScore[1]=0; 

        fread(pucQltyScore,1,1,pfInFile);

        fprintf(pfIn,"%u",pucQltyScore[0]); //fprintf(stdout,"QltyScore-%u",pucQltyScore[0]);
          
        //pucQltyScore will changed to  multibyte in future...
        //for (int i=1; i< ucCnt; i++) fprintf(pfIn,",%u",pucQltyScore[i]);

        fread(&usPEnd,2,1,pfInFile); fprintf(pfIn,"\t%u\t",usPEnd);   //fprintf(stdout,"\tPEnd-%u\t",usPEnd);
        fread(&cStrand,1,1,pfInFile); fprintf(pfIn,"%d\n",cStrand);   //fprintf(stdout,"Strand-%d\n",cStrand);
    }
    fclose(pfIn);  fseek(pfInFile,0,SEEK_SET);
}


void CSXReadSNV::PrintSortedMNV(stMNVRecList aMNVRecList[])
{
    fprintf(stdout,"Printing Sorted MNV recs...\n");

    char acFileName[1024]; sprintf(acFileName,"%s/%s_sorted_MNV.txt",m_pcOutputDir,m_pcSampleID);
    FILE *pf=fopen(acFileName,"w"); stMNVRecList *pMNVRecList; stMNVRec *pMNVRec=NULL;
    unsigned long ulTest1=0, ulTest2=0;

    for (unsigned long i=0; i<24; i++)
    {
        pMNVRecList = &aMNVRecList[i];
        for (int j =0; j < pMNVRecList->nUsed; j++)
        {
            if (++ulTest1==100000){fprintf(stdout,"Total Sorted Recs being printed = %lu X 1.0M\n",++ulTest2);ulTest1=0;}

            pMNVRec = &(pMNVRecList->pMNVRec[j]);  

            fprintf(pf,"%lu\t%u\t%s\t%s\t%u\t%u\t%u\t%d\n",i+1,pMNVRec->unOffset,pMNVRec->pcRefBase,pMNVRec->pcVarBase,
                   pMNVRec->ucQryPos,pMNVRec->pucQltyScore[0],pMNVRec->usPEnd,pMNVRec->cStrand);
        }
    }

    fclose(pf);
}


void CSXReadSNV::PrintClusteredMNVRecs(stMNVClusterRec **lpMNVClustered, size_t merge_list_used)
{
    fprintf(stdout,"Print Clustered MNV Recs ...\n");
    char acClusterFileName[1024]; sprintf(acClusterFileName,"%s/%s_clustered_mnv.txt",m_pcOutputDir,m_pcSampleID);
    FILE *pfCluster=fopen(acClusterFileName,"w");

    unsigned long ulTest1=0, ulTest2=0;

    for (unsigned long i=0; i < merge_list_used+1; i++)
    {
         if (++ulTest1==100000) {fprintf(stdout,"Total Clustered Recs printed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}

         fprintf(pfCluster,"%d\t%u\t%s>%s\t%u\t%lu\t%lu\t%u\n",
                 (*lpMNVClustered)[i].nChromosome,(*lpMNVClustered)[i].unOffset,
                 (*lpMNVClustered)[i].pcRefBase,(*lpMNVClustered)[i].pcVarBase,
                 (*lpMNVClustered)[i].unCnt,(*lpMNVClustered)[i].ulTotalQS,
                 (*lpMNVClustered)[i].ulTotalQP,(*lpMNVClustered)[i].ucMaxQS);
    }

    fclose(pfCluster);
}


void CSXReadSNV::GenerateSNVTblNewMNV(FILE *pfInFile, stMNVClusterRec **lpMNVClustered,
                                      size_t &merge_list_alloc, size_t &merge_list_used)
{

    PrintUnSortedMNV(pfInFile); 

    fprintf(stdout,"Loading, converting and splitting...\n");

    int nChro; unsigned char ucCnt; char *pcAllele=NULL;
    unsigned long ulFileOfs,ulTest1=0,ulTest2=0;
    stMNVRecList aMNVRecList[24]; stMNVRecList *pMNVRecList=NULL;  stMNVRec *pstMNVRec=NULL;   
 
    while (!(feof(pfInFile)))
    {
           nChro=0; fread(&nChro,4,1,pfInFile); //fprintf(stdout,"Chro=%d\t",nChro);

           if (nChro==0) break;           
 
           pMNVRecList = &aMNVRecList[nChro-1];
           
           if (pMNVRecList->nUsed == pMNVRecList->nTotal){
               pMNVRecList->nTotal += MNVREC_SIZE;
               pMNVRecList->pMNVRec = (stMNVRec*)realloc(pMNVRecList->pMNVRec, sizeof(stMNVRec)*pMNVRecList->nTotal);
           }

           pstMNVRec = &(pMNVRecList->pMNVRec[pMNVRecList->nUsed]); pMNVRecList->nUsed++;
                  
           fread(&pstMNVRec->unOffset,4,1,pfInFile); //fprintf(stdout,"Offset=%u\t",pstMNVRec->unOffset);
           fread(&ucCnt,sizeof(ucCnt),1,pfInFile); ulFileOfs = ftell(pfInFile); 

           pcAllele=(char*)malloc(ucCnt+1);

           fread(pcAllele,ucCnt+1,1,pfInFile); 

           if (pcAllele[0]=='-') { //INS
               pstMNVRec->pcRefBase = (char*)malloc(2);   
               pstMNVRec->pcRefBase[0]='-'; pstMNVRec->pcRefBase[1]=0; 
     
               pstMNVRec->pcVarBase = (char*)malloc(ucCnt+1);
               memcpy(pstMNVRec->pcVarBase,pcAllele+1,ucCnt); pstMNVRec->pcVarBase[ucCnt]=0; 

           } else if (pcAllele[ucCnt]=='-'){ //DEL

               pstMNVRec->pcRefBase = (char*)malloc(ucCnt+1);
               memcpy(pstMNVRec->pcRefBase,pcAllele,ucCnt); pstMNVRec->pcRefBase[ucCnt]=0;         

               pstMNVRec->pcVarBase = (char*)malloc(2); 
               pstMNVRec->pcVarBase[0]='-'; pstMNVRec->pcVarBase[1]=0;

           } else { //SNP
               pstMNVRec->pcRefBase = (char*)malloc(ucCnt+1);
               pstMNVRec->pcVarBase = (char*)malloc(ucCnt+1);

               fseek(pfInFile,ulFileOfs,SEEK_SET);
               fread(pstMNVRec->pcRefBase,ucCnt,1,pfInFile); pstMNVRec->pcRefBase[ucCnt]=0;
               fread(pstMNVRec->pcVarBase,ucCnt,1,pfInFile); pstMNVRec->pcVarBase[ucCnt]=0;
           }    
            
           //fprintf(stdout,"RefBase=%s VarBase=%s\t",pstMNVRec->pcRefBase,pstMNVRec->pcVarBase);  

           if (pcAllele) free(pcAllele);          
            
           fread(&(pstMNVRec->ucQryPos),1,1,pfInFile); //fprintf(stdout,"QryPos = %u\t",pstMNVRec->ucQryPos);

           //if (pstMNVRec->pcVarBase[0]=='-'){ //DEL   //pstMNVRec->pucQltyScore will changed to multibyte i future.         
               pstMNVRec->pucQltyScore = (unsigned char*)malloc(2);
               fread(pstMNVRec->pucQltyScore,1,1,pfInFile); pstMNVRec->pucQltyScore[1]=0;               
           /*} else {    
               pstMNVRec->pucQltyScore = (unsigned char*)malloc(ucCnt+1); 
               fread(pstMNVRec->pucQltyScore,ucCnt,1,pfInFile); pstMNVRec->pucQltyScore[ucCnt]=0;
           }*/

           //fprintf(stdout,"QltyScore=%u\t",pstMNVRec->pucQltyScore[0]);  
           fread(&(pstMNVRec->usPEnd),2,1,pfInFile);  //fprintf(stdout,"PEnd=%u\t",pstMNVRec->usPEnd);
           fread(&(pstMNVRec->cStrand),1,1,pfInFile); //fprintf(stdout,"Strand=%d\n",pstMNVRec->cStrand);
    }

    fprintf(stdout,"sorting...\n");

    for (int i=0; i<24; i++){
         if (pthread_create(&(aMNVRecList[i].ThreadID),NULL,SortMNVData,(void*)(&aMNVRecList[i])) != 0)
         { fprintf(stdout,"Failed to spawn SortMNVData Thread %d ...\n",i); exit(0); }
    }

    for (int i=0; i<24; i++){
         if (pthread_join(aMNVRecList[i].ThreadID,NULL)!=0)
         { fprintf(stdout,"Failed to perform Thread Join %d ...\n",i); exit(0); }
    }

    PrintSortedMNV(aMNVRecList);

    fprintf(stdout,"alloc memory for clustering...\n");
    AllocMNVMergeListSize(lpMNVClustered,merge_list_alloc,merge_list_used);

    fprintf(stdout,"clustering...\n");

    char *pcRefBase_Curr=NULL,*pcVarBase_Curr=NULL,cStrand;  
    int nChromosome_Curr;
    unsigned int unOffset_Curr=0;
    unsigned short usQP, usQS, usPEnd=0,ausPEType[4]={49152,32768,16384,0}; 
     
    int nChroIdx;
    for (nChroIdx=0; nChroIdx<24; nChroIdx++){
         if (aMNVRecList[nChroIdx].nUsed < 1) continue;
         pMNVRecList=&aMNVRecList[nChroIdx]; break;
    }

    for (int idx=0; idx < 4; idx++)
    {
       if ((pMNVRecList->pMNVRec[nChroIdx].usPEnd & ausPEType[idx])==ausPEType[idx])
       { usPEnd = pMNVRecList->pMNVRec[nChroIdx].usPEnd - ausPEType[idx]; break; }
    }

    usQP = pMNVRecList->pMNVRec[0].ucQryPos;

    if (pMNVRecList->pMNVRec[0].pucQltyScore[0] < m_nQScoreOfs) usQS = 0;  // Will be invalid if pcQltyScore changed to multibyte
    else usQS = pMNVRecList->pMNVRec[0].pucQltyScore[0]-m_nQScoreOfs;

    nChromosome_Curr = nChroIdx+1;
    cStrand = pMNVRecList->pMNVRec[0].cStrand & 0x01;
    unOffset_Curr = pMNVRecList->pMNVRec[0].unOffset; 

    (*lpMNVClustered)[merge_list_used].nChromosome = nChromosome_Curr;
    (*lpMNVClustered)[merge_list_used].unOffset = unOffset_Curr;

    (*lpMNVClustered)[merge_list_used].pcRefBase = new char[strlen(pMNVRecList->pMNVRec[0].pcRefBase)+1];
    strcpy((*lpMNVClustered)[merge_list_used].pcRefBase,pMNVRecList->pMNVRec[0].pcRefBase);
    pcRefBase_Curr=(*lpMNVClustered)[merge_list_used].pcRefBase;   

    (*lpMNVClustered)[merge_list_used].pcVarBase = new char[strlen(pMNVRecList->pMNVRec[0].pcVarBase)+1];
    strcpy((*lpMNVClustered)[merge_list_used].pcVarBase,pMNVRecList->pMNVRec[0].pcVarBase);
    pcVarBase_Curr=(*lpMNVClustered)[merge_list_used].pcVarBase;

    (*lpMNVClustered)[merge_list_used].ulTotalQP = usQP;
    (*lpMNVClustered)[merge_list_used].ulTotalQS = usQS;
    (*lpMNVClustered)[merge_list_used].ucMaxQS = usQS;
    (*lpMNVClustered)[merge_list_used].unCnt = 1;

    if (usPEnd != 0) (*lpMNVClustered)[merge_list_used].unPEndCnt = 1;

    if (cStrand==0) (*lpMNVClustered)[merge_list_used].unRvsCnt = 1;
    else if (cStrand==1) (*lpMNVClustered)[merge_list_used].unFwdCnt = 1; 

    for (int i=0; i < 24; i++)
    {
        pMNVRecList=&aMNVRecList[i]; 
         
        for (int j=0; j < pMNVRecList->nUsed; j++)
        {
            if (i==nChroIdx && j==0) continue;
            
            usQP = pMNVRecList->pMNVRec[j].ucQryPos;
            if (pMNVRecList->pMNVRec[j].pucQltyScore[0] < m_nQScoreOfs) usQS=0;  // Will be invalid if pcQltyScore changed to multibyte
            else usQS = pMNVRecList->pMNVRec[j].pucQltyScore[0]-m_nQScoreOfs; 

            for (int idx=0; idx < 4; idx++)
            {
                if ((pMNVRecList->pMNVRec[j].usPEnd & ausPEType[idx])==ausPEType[idx])
                { usPEnd = pMNVRecList->pMNVRec[j].usPEnd - ausPEType[idx]; break; }
            }  
 
            cStrand = pMNVRecList->pMNVRec[j].cStrand & 0x01;
  
            if (nChromosome_Curr != i+1||unOffset_Curr != pMNVRecList->pMNVRec[j].unOffset||
                strcmp(pcRefBase_Curr,pMNVRecList->pMNVRec[j].pcRefBase)!= 0||
                strcmp(pcVarBase_Curr,pMNVRecList->pMNVRec[j].pcVarBase) != 0) 
            {
                nChromosome_Curr = i+1; unOffset_Curr = pMNVRecList->pMNVRec[j].unOffset;

                merge_list_used++;

                if (merge_list_used == merge_list_alloc)
                    AllocMNVMergeListSize(lpMNVClustered,merge_list_alloc,merge_list_used);          

                (*lpMNVClustered)[merge_list_used].nChromosome = nChromosome_Curr;
                (*lpMNVClustered)[merge_list_used].unOffset = unOffset_Curr;

                (*lpMNVClustered)[merge_list_used].pcRefBase = new char[strlen(pMNVRecList->pMNVRec[j].pcRefBase)+1];
                strcpy((*lpMNVClustered)[merge_list_used].pcRefBase,pMNVRecList->pMNVRec[j].pcRefBase);
                pcRefBase_Curr=(*lpMNVClustered)[merge_list_used].pcRefBase;

                (*lpMNVClustered)[merge_list_used].pcVarBase = new char[strlen(pMNVRecList->pMNVRec[j].pcVarBase)+1];
                strcpy((*lpMNVClustered)[merge_list_used].pcVarBase,pMNVRecList->pMNVRec[j].pcVarBase);
                pcVarBase_Curr=(*lpMNVClustered)[merge_list_used].pcVarBase;

                (*lpMNVClustered)[merge_list_used].ulTotalQP = usQP;
                (*lpMNVClustered)[merge_list_used].ulTotalQS = usQS;
                (*lpMNVClustered)[merge_list_used].ucMaxQS = usQS;
                (*lpMNVClustered)[merge_list_used].unCnt = 1;
                if (usPEnd != 0) (*lpMNVClustered)[merge_list_used].unPEndCnt = 1;
                if (cStrand==0) (*lpMNVClustered)[merge_list_used].unRvsCnt = 1;
                else if (cStrand==1) (*lpMNVClustered)[merge_list_used].unFwdCnt = 1;                     
            }   
            else
            {
                (*lpMNVClustered)[merge_list_used].ulTotalQP += usQP;
                (*lpMNVClustered)[merge_list_used].ulTotalQS += usQS;

                if ((*lpMNVClustered)[merge_list_used].ucMaxQS < usQS)
                    (*lpMNVClustered)[merge_list_used].ucMaxQS = usQS; 
                
                (*lpMNVClustered)[merge_list_used].unCnt++;

                if (usPEnd != 0) (*lpMNVClustered)[merge_list_used].unPEndCnt++;

                if (cStrand==0) (*lpMNVClustered)[merge_list_used].unRvsCnt++;
                else if (cStrand==1) (*lpMNVClustered)[merge_list_used].unFwdCnt++;   
            } 

            if (++ulTest1 == 1000000){fprintf(stdout,"Total Recs gone through clustering = %luM\n", ++ulTest2); ulTest1=0;}  
        } //end for 
    } //end for

    fprintf(stdout,"Total Clustered Recs = %lu\n", merge_list_used+1);

    if (m_bOutClusteredFile) PrintClusteredMNVRecs(lpMNVClustered,merge_list_used);      
}


void CSXReadSNV::PrintMNVAvgQScoreList()
{
    char acPDensFName[1024]; sprintf(acPDensFName,"%s/PDens_del_%s_mnv.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenMDEL = fopen(acPDensFName,"w"); 

    sprintf(acPDensFName,"%s/PDens_ins_%s_mnv.lst",m_pcOutputDir,m_pcSampleID);
    FILE *pfPDenMINS = fopen(acPDensFName,"w");

    stSNVTblStats astMNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astMINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astMDELStats[m_unMaxSupportingReads+1];

    int nMNPIdx,nMINSIdx,nMDELIdx,nBps;
    unsigned int unChroIdx, unStart, unStop, unReadDensity=0, unPDensity=0;
    unsigned long ulTest1=0, ulTest2=0;
    bool bdbSNP=false; char acRsid[4096];
    float fAvgQryPos,fAvgQScore;

    FILE *pfMNP,*pfMINS,*pfMDEL,*pfMNPStats,*pfMINSStats,*pfMDELStats;
    pfMNP=pfMINS=pfMDEL=pfMNPStats=0,pfMINSStats=pfMDELStats=0;

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acFileName[24+strlen(m_pcSampleID)];
    sprintf(acFileName,"%s/%s_mnp_c1_aqs_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMNP = fopen(acFileName,"w");
    if (!pfMNP) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}


    sprintf(acFileName,"%s/%s_mnp_summary_c1_aqs_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMNPStats = fopen(acFileName,"w");
    if (!pfMNPStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    sprintf(acFileName,"%s/%s_ins_c1_aqs_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMINS = fopen(acFileName,"w");
    if (!pfMINS) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_ins_summary_c1_aqs_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMINSStats = fopen(acFileName,"w");
    if (!pfMINSStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////
 
    sprintf(acFileName,"%s/%s_del_c1_aqs_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMDEL = fopen(acFileName,"w");
    if (!pfMDEL) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_del_summary_c1_aqs_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
            ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMDELStats = fopen(acFileName,"w");
    if (!pfMDELStats) {fprintf(stdout,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

///////

    fprintf(pfMNP,"#REPORT NAME\tMNP Avg. Score Distribution\n");
    fprintf(pfMNP,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMNP,"#LANE NO\tALL\n");
    fprintf(pfMNP,"#GENERATED AT\t"); PrintRptDateTime(pfMNP);
    fprintf(pfMNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfMNP,"#FILTER\t>=1\n\n\n");

    fprintf(pfMNP,"Chromosome\tStart\tStop\tNucleotide_Variant\tFwd_MNP_Reads"
                  "\tRvs_MNP_Reads\tTotal_MNP_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");

    fprintf(pfMINS,"#REPORT NAME\tINS(MNV) Avg. Score Distribution\n");
    fprintf(pfMINS,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMINS,"#LANE NO\tALL\n");
    fprintf(pfMINS,"#GENERATED AT\t"); PrintRptDateTime(pfMINS);
    fprintf(pfMINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfMINS,"#FILTER\t>=1\n\n\n");

    fprintf(pfMINS,"Chromosome\tStart\tStop\tInserted_Base\tFwd_INS(MNV)_Reads"
                   "\tRvs_INS(MNV)_Reads\tTotal_INS(MNV)_Reads"
                   "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");

    fprintf(pfMDEL,"#REPORT NAME\tDEL(MNV) Avg. Score Distribution\n");
    fprintf(pfMDEL,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMDEL,"#LANE NO\tALL\n");
    fprintf(pfMDEL,"#GENERATED AT\t"); PrintRptDateTime(pfMDEL);
    fprintf(pfMDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);
    fprintf(pfMDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfMDEL,"Chromosome\tStart\tStop\tDeleted_Base\tFwd_DEL(MNV)_Reads"
                  "\tRvs_DEL(MNV)_Reads\tTotal_DEL(MNV)_Reads"
                  "\tRead_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP\n");
       
    for (unsigned long i=0; i < (m_merge_list_used+1); i++)
    {
         if (++ulTest1 == 1000000){fprintf(stdout,"Total Clustered Recs Processed = %luM\n", ++ulTest2); ulTest1=0; }

         if (m_pMNVClustered[i].nChromosome == 0) continue;

         unChroIdx = m_pMNVClustered[i].nChromosome-1;

         if (m_pMNVClustered[i].pcRefBase[0]!='-' && m_pMNVClustered[i].pcVarBase[0]!='-') 
         {
             bdbSNP=OutputMNPNovel(m_pMNVClustered[i].nChromosome,m_pMNVClustered[i].unOffset,
                                   m_pMNVClustered[i].pcVarBase,pfMNP,acRsid);

             fAvgQryPos = CalcAvgQS(m_pMNVClustered[i].ulTotalQP,m_pMNVClustered[i].unCnt);
             fAvgQScore = CalcAvgQS(m_pMNVClustered[i].ulTotalQS,m_pMNVClustered[i].unCnt);             

             nBps = strlen(m_pMNVClustered[i].pcRefBase); 
             unStart = m_pMNVClustered[i].unOffset;
             unStop = (unStart+nBps)-1;

             for (unsigned int unOffset=unStart; unOffset < (unStop+1); unOffset++)
             { unReadDensity += GetReadDensity(m_apfRD[unChroIdx], unOffset); }
            
             unReadDensity =  (unsigned int)(((double)unReadDensity/(double)nBps)+0.5) + m_pMNVClustered[i].unCnt; 

             if (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
             if (Percentage(m_pMNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

             OutputChrom(pfMNP,m_pMNVClustered[i].nChromosome);

             fprintf(pfMNP,"\t%u\t%u\t%s>%s\t%u\t%u\t%u",unStart,unStop,
                     m_pMNVClustered[i].pcRefBase,m_pMNVClustered[i].pcVarBase,m_pMNVClustered[i].unFwdCnt,
                     m_pMNVClustered[i].unRvsCnt,m_pMNVClustered[i].unCnt);
            
             fprintf(pfMNP,"\t%u\t%u\t%.2f\t%.2f\t%s",
                     unReadDensity,m_pMNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

             fprintf(pfMNP,"\n");

             nMNPIdx = (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pMNVClustered[i].unCnt-1;
             astMNPStats[nMNPIdx].unCnt++;
             if (bdbSNP) astMNPStats[nMNPIdx].undbSNP++;
             astMNPStats[nMNPIdx].ulTotalReadDens += unReadDensity; 
         }   
         else if (m_pMNVClustered[i].pcRefBase[0] =='-' )
         {
             bdbSNP=m_pObjdbINS->OutputMNVINDELNovel(m_pMNVClustered[i].nChromosome,m_pMNVClustered[i].unOffset,
                                                     m_pMNVClustered[i].pcVarBase, pfMINS,acRsid);   

             fAvgQryPos = CalcAvgQS(m_pMNVClustered[i].ulTotalQP,m_pMNVClustered[i].unCnt);
             fAvgQScore = CalcAvgQS(m_pMNVClustered[i].ulTotalQS,m_pMNVClustered[i].unCnt);  

             nBps = strlen(m_pMNVClustered[i].pcVarBase);
             unStart = m_pMNVClustered[i].unOffset;
             unStop = (unStart+nBps)-1;
             unReadDensity = 0;  

             for (unsigned int unOffset=unStart; unOffset < (unStop+1); unOffset++)               
             { unReadDensity += GetReadDensityEx(m_apfRD[unChroIdx], unOffset); }
             unReadDensity = (unsigned int)(((double)unReadDensity/(double)nBps)+0.5);

             if (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
             if (Percentage(m_pMNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

             OutputChrom(pfMINS,m_pMNVClustered[i].nChromosome);

             fprintf(pfMINS,"\t%u\t%u\t%s\t%u\t%u\t%u",unStart,unStop,  
                     m_pMNVClustered[i].pcVarBase,m_pMNVClustered[i].unFwdCnt,
                     m_pMNVClustered[i].unRvsCnt,m_pMNVClustered[i].unCnt);

             fprintf(pfMINS,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                     unReadDensity,m_pMNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

             fprintf(pfMINS,"\n");

             nMINSIdx = (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pMNVClustered[i].unCnt-1;
             astMINSStats[nMINSIdx].unCnt++;
             if (bdbSNP) astMINSStats[nMINSIdx].undbSNP++;
             astMINSStats[nMINSIdx].ulTotalReadDens += unReadDensity;

             unPDensity = 0; 
             for (unsigned int unOffset=unStart; unOffset < (unStop+1); unOffset++)
             { unPDensity = GetReadDensityEx(m_apfPDens[unChroIdx], unOffset); }
             unPDensity = (unsigned int)(((double)unPDensity/(double)nBps)+0.5);

             OutputChrom(pfPDenMINS,m_pMNVClustered[i].nChromosome);

             fprintf(pfPDenMINS,"\t%u\t%u\t%u\n",unStart,unStop,unPDensity+unReadDensity);
         } 
         else
         {
             bdbSNP=m_pObjdbDEL->OutputMNVINDELNovel(m_pMNVClustered[i].nChromosome,m_pMNVClustered[i].unOffset,
                                                     m_pMNVClustered[i].pcRefBase,pfMDEL,acRsid);

             fAvgQryPos = CalcAvgQS(m_pMNVClustered[i].ulTotalQP,m_pMNVClustered[i].unCnt);
             fAvgQScore = CalcAvgQS(m_pMNVClustered[i].ulTotalQS,m_pMNVClustered[i].unCnt);

             nBps = strlen(m_pMNVClustered[i].pcRefBase);
             unStart = m_pMNVClustered[i].unOffset;
             unStop = (unStart+nBps)-1;
             unReadDensity = 0;

             for (unsigned int unOffset=unStart; unOffset < (unStop+1); unOffset++)
             { unReadDensity += GetReadDensity(m_apfRD[unChroIdx], unOffset); }
             unReadDensity = (unsigned int)(((double)unReadDensity/(double)nBps)+0.5) + m_pMNVClustered[i].unCnt;

             if (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
             if (Percentage(m_pMNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

             OutputChrom(pfMDEL,m_pMNVClustered[i].nChromosome);
             fprintf(pfMDEL,"\t%u\t%u\t%s\t%u\t%u\t%u",unStart,unStop,
                     m_pMNVClustered[i].pcRefBase,m_pMNVClustered[i].unFwdCnt,
                     m_pMNVClustered[i].unRvsCnt,m_pMNVClustered[i].unCnt);

             fprintf(pfMDEL,"\t%u\t%u\t%.2f\t%.2f\t%s\t",
                     unReadDensity,m_pMNVClustered[i].unPEndCnt,fAvgQryPos,fAvgQScore,acRsid);

             fprintf(pfMDEL,"\n");

             nMDELIdx = (m_pMNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pMNVClustered[i].unCnt-1;
             astMDELStats[nMDELIdx].unCnt++;
             if (bdbSNP) astMDELStats[nMDELIdx].undbSNP++;
             astMDELStats[nMDELIdx].ulTotalReadDens += unReadDensity;

             unPDensity = 0;
             for (unsigned int unOffset=unStart; unOffset < (unStop+1); unOffset++)
             { unPDensity = GetReadDensity(m_apfPDens[unChroIdx], unOffset); }
             unPDensity = (unsigned int)(((double)unPDensity/(double)nBps)+0.5) + m_pMNVClustered[i].unCnt;

             OutputChrom(pfPDenMDEL,m_pMNVClustered[i].nChromosome);

             fprintf(pfPDenMDEL,"\t%u\t%u\t%u\n",unStart,unStop,unPDensity+unReadDensity);

         }
    } //end for 

    fprintf(stderr,"Printing AvgQScore Statistic Report...\n");

    PrintAvgQSStatsRpt(pfMNPStats,astMNPStats,(char*)"MNP");
    PrintAvgQSStatsRpt(pfMINSStats,astMINSStats,(char*)"INS(MNV)");
    PrintAvgQSStatsRpt(pfMDELStats,astMDELStats,(char*)"DEL(MNV)");

ExitFunc:
    if (pfMNP) fclose(pfMNP); if (pfMINS) fclose(pfMINS);
    if (pfMDEL) fclose(pfMDEL); if (pfMNPStats) fclose(pfMNPStats);
    if (pfMINSStats) fclose(pfMINSStats); if (pfMDELStats) fclose(pfMDELStats);
    if (pfPDenMINS) fclose(pfPDenMINS); if (pfPDenMDEL) fclose(pfPDenMDEL);
}


void CSXReadSNV::OutputAvgQScoreTableMNV()
{
    FILE *pfMNV=NULL, *pfdbSnp=NULL, *pfdbIndel=NULL, *pfdbINS=NULL, *pfdbDEL=NULL;

    char acFile[strlen(m_pcDensFilePath)+16];
    m_apfRD = new FILE*[24]; m_apfPDens = new FILE*[24];
    for (int i=0; i<24; i++){ m_apfRD[i]=NULL; m_apfPDens[i]=NULL; } 

    pfMNV = fopen(m_pcInFile,"rb"); if (!pfMNV){fprintf(stdout,"Failed to open %s ...\n",m_pcInFile); goto Exit;} 

    pfdbSnp = fopen(m_pcdbSnp,"r"); if (!pfdbSnp){fprintf(stdout,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;}

    if (m_pcdbIndel){
        pfdbIndel = fopen(m_pcdbIndel,"r"); if (!pfdbSnp){fprintf(stdout,"Failed to open %s ...\n",m_pcdbSnp); goto Exit;}
    }
    else {
       pfdbINS = fopen(m_pcdbINS,"r"); if (!pfdbINS) {fprintf(stderr,"Failed to open %s ...\n",m_pcdbINS); goto Exit;}
       pfdbDEL = fopen(m_pcdbDEL,"r"); if (!pfdbDEL) {fprintf(stderr,"Failed to open %s ...\n",m_pcdbDEL); goto Exit;}
    }

    for (int i=0; i<24;i++)
    {
        sprintf(acFile,"%s/%s.seq%d.read_den_all",m_pcDensFilePath,m_pcProjSampleID,i+1);
        m_apfRD[i] = fopen(acFile,"r");
        if (!m_apfRD[i]){fprintf(stdout,"Failed to open %s ...\n", acFile); goto Exit;}
        else fprintf(stdout,"Touch %s ...\n",acFile);
    }

    for (int i=0; i<24;i++)
    {
        sprintf(acFile,"%s/%s.seq%d.read_den_perfect",m_pcDensFilePath,m_pcProjSampleID,i+1);
        m_apfPDens[i] = fopen(acFile,"r");
        if (!m_apfPDens[i]){fprintf(stdout,"Failed to open %s ...\n", acFile); goto Exit;}
        else fprintf(stdout,"Touch %s ...\n",acFile);
    }

    fprintf(stdout,"Generating %s Map...\n",m_pcdbSnp);

    if (!GenerateDBSNPMap(pfdbSnp)) goto Exit;

    fprintf(stdout,"Generating %s Map...\n",m_pcdbINS);
    m_pObjdbINS = new CSXDBIndelgChecker(); if (!m_pObjdbINS) goto Exit;
    if (!m_pObjdbINS->GenerateDBINDELMap(pfdbINS)) goto Exit;

    fprintf(stdout,"Generating %s Map...\n",m_pcdbDEL);
    m_pObjdbDEL = new CSXDBIndelgChecker(); if (!m_pObjdbDEL) goto Exit;
    if (!m_pObjdbDEL->GenerateDBINDELMap(pfdbDEL)) goto Exit;

    fprintf(stdout,"Generating MNV Table...\n");

    GenerateSNVTblNewMNV(pfMNV,&m_pMNVClustered,m_merge_list_alloc,m_merge_list_used);

    fprintf(stdout,"Printing MNP, INS(MNV), DEL(MNV) Lists...\n"); PrintMNVAvgQScoreList();

Exit:
    if (pfMNV) fclose(pfMNV);
    if (pfdbSnp) fclose(pfdbSnp); if (pfdbIndel) fclose(pfdbIndel);
    ClrDBSNPMap(); ClrDBINDELMap();  

    for (int j=0;j<24; j++){
         if (m_apfRD[j]) fclose(m_apfRD[j]);
         if (m_apfPDens[j]) fclose(m_apfPDens[j]); 
    }

    delete[] m_apfRD; delete[] m_apfPDens;
}


void CSXReadSNV::OutputMNVList()
{
    FILE *pfSrc=NULL;
    
    pfSrc = fopen(m_pcInFile,"r"); if (!pfSrc) {fprintf(stdout,"Failed to open %s ...\n",m_pcInFile); goto Exit;}
    fprintf(stdout,"Loading Annotation Data Source...\n");

    m_pSXAnnotate = new SXAnnotate;
    
    if (!m_pSXAnnotate->LoadDataSource(m_pcBinnedFile,m_pcExceptionFile,m_pcAnnotateFile)) goto Exit;

    m_punGINums = strcmp(m_pcGenomeVer,"37.1")==0?g_unGINums_G37:g_unGINums_G36;
    
    if (m_eMNV==eMNP){
        fprintf(stdout,"Printing MNP List Table...\n");PrintMNPList(pfSrc);
    }
    else if (m_eMNV==eMINS){
        fprintf(stdout,"Printing INS(MNV) List Table...\n");PrintMINSList(pfSrc);
    }
    else if (m_eMNV==eMDEL){
        fprintf(stdout,"Printing DEL(MNV) List Table...\n");PrintMDELList(pfSrc);
    }
    
Exit:

    if (pfSrc) fclose(pfSrc);
    if (m_pSXAnnotate) {delete m_pSXAnnotate; m_pSXAnnotate=NULL;}
}


void CSXReadSNV::PrintMNPList(FILE* pfSrc)
{
    char acbuf[4096],*pcAllele=NULL,acdbSnp[1024],acHomoHet[10],cChromosome=0;
    char *pChr=NULL; float fLocalCNV=0.0,fAvgQScore=0.0,fAvgQryPos=0.0;
    unsigned int unChromosome=0,unGINum=0;
    unsigned long ulStart=0,ulStop=0,ulFwdSR=0,ulRvsSR=0,ulMNPReads=0,ulReadDensity=0,ulPEnds=0;
    unsigned long ulMappableBases=0, ulMappableReadDens=0, ulMappableReptDens=0;
    unsigned long ulMappableReptCnt=0, ulMappableReptVal=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nMNPIdx=0; stSNVTblStats astMNPStats[m_unMaxSupportingReads+1];

    FILE *pfMNP,*pfMNPStats; pfMNP=pfMNPStats=0;
    
    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)],*pRecs=NULL; stAnnotate aAnnotate;    
    //int nTotalTabs=0; 

    sprintf(acFileName,"%s/%s_mnp_c1_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMNP = fopen(acFileName,"w");
    if (!pfMNP) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_mnp_summary_c1_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMNPStats = fopen(acFileName,"w");
    if (!pfMNPStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pfMNP,"#REPORT NAME\tMNP\n");
    fprintf(pfMNP,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMNP,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMNP,"#LANE NO\tALL\n");
    fprintf(pfMNP,"#GENERATED AT\t"); PrintRptDateTime(pfMNP);
    fprintf(pfMNP,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMNP,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);

    fprintf(pfMNP,"\n\n\n");

    fprintf(pfMNP,"Chromosome\tGiNumber\tStart\tStop\tNucleotide_Variant\tFwd_MNP_Reads"
                  "\tRvs_MNP_Reads\tTotal_MNP_Reads\tTotal_Read_Density\tPEnd_Count"
                  "\tAvg_QryPos\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                  "\tZygosity\tLocal_Copy_Number\tMappable_Bases"
                  "\tMappable_Read_Densities\tMappable_Repeat_Densities"
                  "\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                  "\tKnown_CNVRegion\n");

    while (!feof(pfSrc))
    {
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' || acbuf[0]==0 || acbuf[0]=='C') continue;
/*
        if (nTotalTabs == 0){
           pChr = strchr(acbuf,'\t');
           while(pChr){ nTotalTabs++; pChr = strchr(pChr+1,'\t');}
        }
*/
        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfMNP,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if(pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        pChr = strtok(NULL,"\t"); ulStart = atol(pChr);           //Start
        pChr = strtok(NULL,"\t"); ulStop = atol(pChr);            //Stop

        pChr = strtok(NULL,"\t");                                 //Nucleotide_Variant
        pcAllele=new char[strlen(pChr)+1]; strcpy(pcAllele,pChr); //Nucleotide_Variant

        pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);             //Fwd_MNP_Reads
        pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);             //Rvs_MNP_Reads 
        pChr = strtok(NULL,"\t"); ulMNPReads=atol(pChr);          //MNP_Reads
        pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr);       //Read_Density
        pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);             //PEnd_Count
        pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);         //Avg_QryPos
        pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);          //Avg_QScore
        pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);           //dbSnp
        pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);         //Homo_Het
        pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr);         //Local_CNV
        pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr);     //Mappable_Bases
        pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr);  //Mappable_Read_Densies
        pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr);  //Mappale_Repeat_Densities
        pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr);   //Mappable_Bases_Repeat
        pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr);   //Repeat_Densities_MBP

        unGINum = m_punGINums[unChromosome-1];
        fprintf(pfMNP,"\t%u\t%lu\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulStart,ulStop,pcAllele,ulFwdSR,ulRvsSR,ulMNPReads,ulReadDensity);

        fprintf(pfMNP,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

        pRecs = m_pSXAnnotate->GetRecord(cChromosome,ulStart,ulStop);

        ExtractAnnotation(pRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfMNP);fprintf(pfMNP,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfMNP);fprintf(pfMNP,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfMNP); fprintf(pfMNP,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfMNP); fprintf(pfMNP,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfMNP,"Y\t");
        else fprintf(pfMNP,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfMNP,"Y\t");
        else fprintf(pfMNP,"N\t");

        //PrintVECListContents(aAnnotate.UTR_List, pfMNP); fprintf(pfMNP,"\t");  

        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfMNP);fprintf(pfMNP,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfMNP,"-\t");}
        
        fprintf(pfMNP,"%s\t",acHomoHet);        
        fprintf(pfMNP,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
                fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                ulMappableReptCnt,ulMappableReptVal);
        
        if (aAnnotate.CNV_List.size() > 0) fprintf(pfMNP,"Y");
        else fprintf(pfMNP,"N");

        fprintf(pfMNP,"\n");

        nMNPIdx = (ulMNPReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulMNPReads-1;

        astMNPStats[nMNPIdx].unCnt++;
        astMNPStats[nMNPIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astMNPStats[nMNPIdx].undbSNP++;
        if (bGene) astMNPStats[nMNPIdx].unGene++;
        if (bExon) astMNPStats[nMNPIdx].unExon++;       

        aAnnotate.ClrList(); if(pcAllele) delete[] pcAllele;         

    } //end while

    fprintf(stdout,"Printing MNP Statistic Report...\n");

    if (pfMNPStats) PrintStatsRpt(pfMNPStats,astMNPStats,(char*)"MNP");

ExitFunc:
    if (pfMNP) fclose(pfMNP); if (pfMNPStats) fclose(pfMNPStats);
}


void CSXReadSNV::PrintMINSList(FILE* pfSrc)
{
    char acbuf[4096],*pcAllele=NULL,acdbSnp[1024],acHomoHet[10]; char *pChr=NULL, cChromosome=0;
    float fLocalCNV=0.0,fAvgQScore=0.0,fAvgQryPos=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulStart=0,ulStop=0,ulMINSReads=0,ulFwdSR=0,ulRvsSR=0,ulReadDensity=0,ulPEnds=0;
    unsigned long ulMappableBases=0,ulMappableReadDens=0,ulMappableReptDens=0;
    unsigned long ulMappableReptCnt=0,ulMappableReptVal=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nMINSIdx; stSNVTblStats astMINSStats[m_unMaxSupportingReads+1];

    FILE *pfMINS,*pfMINSStats; pfMINS=pfMINSStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)],*pcRecs=NULL; stAnnotate aAnnotate; 
    //int nTotalTabs=0;

    sprintf(acFileName,"%s/%s_ins_c1_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMINS = fopen(acFileName,"w");
    if (!pfMINS) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_ins_summary_c1_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMINSStats = fopen(acFileName,"w");
    if (!pfMINSStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    m_stSNVAnalysis.Init();

    fprintf(pfMINS,"#REPORT NAME\tInsertion(MNV)\n");
    fprintf(pfMINS,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMINS,"#LANE NO\tALL\n");
    fprintf(pfMINS,"#GENERATED AT\t"); PrintRptDateTime(pfMINS);
    fprintf(pfMINS,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMINS,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);

    fprintf(pfMINS,"\n\n\n");

    fprintf(pfMINS,"Chromosome\tGiNumber\tStart\tStop\tInserted_Base\tFwd_INS(MNV)_Reads\tRvs_INS(MNV)_Reads"
                   "\tTotal_INS(MNV)_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                   "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                   "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                   "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                   "\tKnown_CNV_Region\n");

    while (!feof(pfSrc))
    {
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2); ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;} pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#' ||acbuf[0]==0|| acbuf[0]=='C') continue;
        /* 
        if (nTotalTabs == 0) {
           pChr = strchr(acbuf,'\t');
           while(pChr){ nTotalTabs++; pChr = strchr(pChr+1,'\t');}
        }
        */
        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfMINS,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if (pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        pChr = strtok(NULL,"\t"); ulStart = atol(pChr);               //Start
        pChr = strtok(NULL,"\t"); ulStop = atol(pChr);                //Stop

        pChr = strtok(NULL,"\t");                                     //Inserted base
        pcAllele = new char[strlen(pChr)+1]; strcpy(pcAllele,pChr);

        pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);                 //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);                 //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t"); ulMINSReads=atol(pChr);             //SNP_Reads
        pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr);           //Read_Density
        pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);                 //PEnd_Count
        pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);             //Avg_QryPos
        pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);              //Avg_QScore
        pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);               //dbSnp
        pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);             //Homo_Het
        pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr);             //Local CNV
        pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr);         //Mappable_Bases
        pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr);      //Mappable_Read_Densities
        pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr);      //Mappable_Repeat_Densities
        pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr);       //Mappable_Bases_Repeat
        pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr);       //Repeat_Densities_MBP

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfMINS,"\t%u\t%lu\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulStart,ulStop,pcAllele,ulFwdSR,ulRvsSR,ulMINSReads,ulReadDensity);

        fprintf(pfMINS,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

        pcRecs = m_pSXAnnotate->GetRecord(cChromosome,ulStart,ulStop);
        ExtractAnnotation(pcRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfMINS);fprintf(pfMINS,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfMINS);fprintf(pfMINS,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfMINS); fprintf(pfMINS,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfMINS); fprintf(pfMINS,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfMINS,"Y\t");
        else fprintf(pfMINS,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfMINS,"Y\t");
        else fprintf(pfMINS,"N\t");
 
        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfMINS);fprintf(pfMINS,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfMINS,"-\t");}

        fprintf(pfMINS,"%s\t",acHomoHet);
        fprintf(pfMINS,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
                   fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                   ulMappableReptCnt,ulMappableReptVal);

        if (aAnnotate.CNV_List.size() > 0) fprintf(pfMINS,"Y\n");
        else fprintf(pfMINS,"N\n");

        nMINSIdx = (ulMINSReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulMINSReads-1;

        astMINSStats[nMINSIdx].unCnt++;
        astMINSStats[nMINSIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astMINSStats[nMINSIdx].undbSNP++;
        if (bGene) astMINSStats[nMINSIdx].unGene++;
        if (bExon) astMINSStats[nMINSIdx].unExon++;
        aAnnotate.ClrList();

        /*  
        if (strcmp(acHomoHet,"het")==0) m_stSNVAnalysis.ulTotalSNPHet++;
        else if (strcmp(acHomoHet,"hom")==0) m_stSNVAnalysis.ulTotalSNPHom++;
        else m_stSNVAnalysis.ulTotalSNPGreyArea++;
        */ 
        if(pcAllele) delete[] pcAllele;
   } //end while

   fprintf(stdout,"Printing MINS Statistic Report...\n");

   if (pfMINSStats) PrintStatsRpt(pfMINSStats,astMINSStats,(char*)"MINS");

ExitFunc:
   if (pfMINS) fclose(pfMINS); if (pfMINSStats) fclose(pfMINSStats);
}


void CSXReadSNV::PrintMDELList(FILE* pfSrc)
{
    char acbuf[4096],*pcAllele,acdbSnp[1024]; char *pChr=NULL,acHomoHet[10],cChromosome;
    float fLocalCNV=0.0,fAvgQryPos=0.0,fAvgQScore=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulStart=0,ulStop=0,ulFwdSR=0,ulRvsSR=0,ulMDELReads=0,ulPEnds=0,ulReadDensity=0;
    unsigned long ulMappableBases=0,ulMappableReadDens=0,ulMappableReptDens=0;
    unsigned long ulMappableReptCnt=0,ulMappableReptVal=0,ulTest1=0,ulTest2=0;

    bool bGene,bExon; int nMDELIdx; stSNVTblStats astMDELStats[m_unMaxSupportingReads+1];

    FILE *pfMDEL,*pfMDELStats; pfMDEL=pfMDELStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)],*pcRecs=NULL; stAnnotate aAnnotate; //int nTotalTabs=0;

    sprintf(acFileName,"%s/%s_del_c1_%02d%02d%02d_mnv.lst",m_pcOutputDir,m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMDEL = fopen(acFileName,"w");
    if (!pfMDEL) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    sprintf(acFileName,"%s/%s_del_summary_c1_%02d%02d%02d_mnv.rpt",m_pcOutputDir,m_pcSampleID,
           ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    pfMDELStats = fopen(acFileName,"w");
    if (!pfMDELStats) {fprintf(stderr,"Failed to open %s ...\n",acFileName); goto ExitFunc;}

    fprintf(pfMDEL,"#REPORT NAME\tDeletion(MNV)\n");
    fprintf(pfMDEL,"#PROJECT NAME\t%s\n",m_pcProjSampleID);
    fprintf(pfMDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfMDEL,"#LANE NO\tALL\n");
    fprintf(pfMDEL,"#GENERATED AT\t"); PrintRptDateTime(pfMDEL);
    fprintf(pfMDEL,"#PROGRAM & BUILD\tReadSNV Rev %s\n",VERSION);
    fprintf(pfMDEL,"#REFERENCE DATABASE\tNCBI Human Genome V%s\n",m_pcGenomeVer);

    fprintf(pfMDEL,"\n\n\n");
 
    fprintf(pfMDEL,"Chromosome\tGiNumber\tStart\tStop\tDeleted_Base\tFwd_DEL(MNV)_Reads\tRvs_DEL(MNV)_Reads"
                   "\tTotal_DEL(MNV)_Reads\tTotal_Read_Density\tPEnd_Count\tAvg_QryPos\tAvg_QScore\tdbSNP"
                   "\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tUTR\tExon"
                   "\tZygosity\tLocal_Copy_Number\tMappable_Bases\tMappable_Read_Densities"
                   "\tMappable_Repeat_Densities\tMappable_Bases_Repeat\tRepeat_Densities_MBP"
                   "\tKnown_CNV_Region\n");

    m_stSNVAnalysis.Init();

    while (!feof(pfSrc))
    {
        if (++ulTest1 == 100000){fprintf(stdout,"Total Recs Processed = %lu X 0.1M\n",++ulTest2);ulTest1=0;}

        if (!fgets(acbuf, sizeof(acbuf),pfSrc)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}

        if (acbuf[0]=='#'||acbuf[0]==0||acbuf[0]=='C') continue;
        /* 
        if (nTotalTabs == 0){
           pChr = strchr(acbuf,'\t');
           while(pChr){nTotalTabs++; pChr = strchr(pChr+1,'\t');}
        }
        */
        pChr = strtok(acbuf,"\t"); //Chromosome

        fprintf(pfMDEL,"%s",pChr);

        if (!isalpha(pChr[0])) {unChromosome = atoi(pChr); cChromosome=unChromosome;}
        else if (pChr[0]=='X') {unChromosome = 23; cChromosome='x';}
        else if (pChr[0]=='Y') {unChromosome = 24; cChromosome='y';}
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        pChr = strtok(NULL,"\t"); ulStart = atol(pChr);             //Start
        pChr = strtok(NULL,"\t"); ulStop = atol(pChr);              //Stop 

        pChr = strtok(NULL,"\t");                                   //Deleted_Variant
        pcAllele = new char[strlen(pChr)+1]; strcpy(pcAllele,pChr);
         
        pChr = strtok(NULL,"\t"); ulFwdSR=atol(pChr);               //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t"); ulRvsSR=atol(pChr);               //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t"); ulMDELReads=atol(pChr);           //DEL_Reads
        pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr);         //Read_Density
        pChr = strtok(NULL,"\t"); ulPEnds=atol(pChr);               //PEnd_Count
        pChr = strtok(NULL,"\t"); fAvgQryPos= atof(pChr);           //Avg_QryPos
        pChr = strtok(NULL,"\t"); fAvgQScore=atof(pChr);            //Avg_QScore
        pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr);             //dbSnp
        pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr);           //Homo_Het
        pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr);           //Local CNV
        pChr = strtok(NULL,"\t"); ulMappableBases=atol(pChr);       //Mappable_Bases
        pChr = strtok(NULL,"\t"); ulMappableReadDens=atol(pChr);    //Mappable_Read_Densities
        pChr = strtok(NULL,"\t"); ulMappableReptDens=atol(pChr);    //Mappable_Rept_Densities
        pChr = strtok(NULL,"\t"); ulMappableReptCnt=atol(pChr);     //Mappable_Bases_Repeat
        pChr = strtok(NULL,"\t"); ulMappableReptVal=atol(pChr);     //Repeat_Densities_MBP 

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfMDEL,"\t%u\t%lu\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t",
                unGINum,ulStart,ulStop,pcAllele,ulFwdSR,ulRvsSR,ulMDELReads,ulReadDensity);

        fprintf(pfMDEL,"%lu\t%.2f\t%.2f\t%s\t",ulPEnds,fAvgQryPos,fAvgQScore,acdbSnp);

        pcRecs = m_pSXAnnotate->GetRecord(cChromosome,ulStart,ulStop);
        ExtractAnnotation(pcRecs,aAnnotate);

        PrintVECListContents(aAnnotate.GeneName_List, pfMDEL);fprintf(pfMDEL,"\t");
        bGene = (bool)aAnnotate.GeneName_List.size();

        PrintVECListContents(aAnnotate.GeneDesc_List, pfMDEL);fprintf(pfMDEL,"\t");
        PrintVECListContents(aAnnotate.GeneKW_List, pfMDEL); fprintf(pfMDEL,"\t");
        PrintVECListContents(aAnnotate.miRNA_List, pfMDEL); fprintf(pfMDEL,"\t");

        if (aAnnotate.Promoter_List.size() > 0) fprintf(pfMDEL,"Y\t");
        else fprintf(pfMDEL,"N\t");

        if (aAnnotate.UTR_List.size() > 0) fprintf(pfMDEL,"Y\t");
        else fprintf(pfMDEL,"N\t");

        //PrintVECListContents(aAnnotate.UTR_List, pfMDEL); fprintf(pfMDEL,"\t");

        if (bGene){
            PrintVECListContents(aAnnotate.ExonID_List, pfMDEL);fprintf(pfMDEL,"\t");
            bExon = (bool)aAnnotate.ExonID_List.size();
        }
        else{ bExon = false; fprintf(pfMDEL,"-\t");}

        fprintf(pfMDEL,"%s\t",acHomoHet); 
	fprintf(pfMDEL,"%.2f\t%lu\t%lu\t%lu\t%lu\t%lu\t",
                fLocalCNV,ulMappableBases,ulMappableReadDens,ulMappableReptDens,
                ulMappableReptCnt,ulMappableReptVal);

        if (aAnnotate.CNV_List.size() > 0) fprintf(pfMDEL,"Y\n");
        else fprintf(pfMDEL,"N\n");

        nMDELIdx = (ulMDELReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulMDELReads-1;

        astMDELStats[nMDELIdx].unCnt++;
        astMDELStats[nMDELIdx].ulTotalReadDens += ulReadDensity;
     
        if (acdbSnp[0]!='-') astMDELStats[nMDELIdx].undbSNP++;
        if (bGene) astMDELStats[nMDELIdx].unGene++;
        if (bExon) astMDELStats[nMDELIdx].unExon++;
        aAnnotate.ClrList();

        if (pcAllele) delete[] pcAllele;
         
        //if (strcmp(acHomoHet,"het")==0) m_stSNVAnalysis.ulTotalSNPHet++;
        //else if (strcmp(acHomoHet,"hom")==0) m_stSNVAnalysis.ulTotalSNPHom++;
        //else m_stSNVAnalysis.ulTotalSNPGreyArea++;
         
    } //while

    fprintf(stdout,"Printing DEL(MNV) Statistic Report...\n");

    if (pfMDELStats) PrintStatsRpt(pfMDELStats,astMDELStats,(char*)"DEL(MNV)");

ExitFunc:
    if (pfMDEL) fclose(pfMDEL); if (pfMDELStats) fclose(pfMDELStats);
}

