/* 
 * File:   CSXReadSNV.cc
 * Author: Manish
 * 
 * Created on November 17, 2009, 12:40 PM
 */

#include "CSXReadSNV.h"
#include "SXGeneMapChr.h"
#include "SXExonMapChr.h"
#include "SXDBSNPgChecker.h"
#include "SXDBIndelgChecker.h"
#include "SXCNVMap.h"

#define _TESTING 0

/*
>gi|149288852|ref|NC_000067.5|NC_000067 Mus musculus chromosome 1, reference assembly (C57BL/6J)
>gi|149338249|ref|NC_000068.6|NC_000068 Mus musculus chromosome 2, reference assembly (C57BL/6J)
>gi|149352351|ref|NC_000069.5|NC_000069 Mus musculus chromosome 3, reference assembly (C57BL/6J)
>gi|149354223|ref|NC_000070.5|NC_000070 Mus musculus chromosome 4, reference assembly (C57BL/6J)
>gi|149354224|ref|NC_000071.5|NC_000071 Mus musculus chromosome 5, reference assembly (C57BL/6J)
>gi|149361431|ref|NC_000072.5|NC_000072 Mus musculus chromosome 6, reference assembly (C57BL/6J)
>gi|149361432|ref|NC_000073.5|NC_000073 Mus musculus chromosome 7, reference assembly (C57BL/6J)
>gi|149361523|ref|NC_000074.5|NC_000074 Mus musculus chromosome 8, reference assembly (C57BL/6J)
>gi|149361524|ref|NC_000075.5|NC_000075 Mus musculus chromosome 9, reference assembly (C57BL/6J)
>gi|149288869|ref|NC_000076.5|NC_000076 Mus musculus chromosome 10, reference assembly (C57BL/6J)
>gi|149288871|ref|NC_000077.5|NC_000077 Mus musculus chromosome 11, reference assembly (C57BL/6J)
>gi|149292731|ref|NC_000078.5|NC_000078 Mus musculus chromosome 12, reference assembly (C57BL/6J)
>gi|149292733|ref|NC_000079.5|NC_000079 Mus musculus chromosome 13, reference assembly (C57BL/6J)
>gi|149292735|ref|NC_000080.5|NC_000080 Mus musculus chromosome 14, reference assembly (C57BL/6J)
>gi|149301884|ref|NC_000081.5|NC_000081 Mus musculus chromosome 15, reference assembly (C57BL/6J)
>gi|149304713|ref|NC_000082.5|NC_000082 Mus musculus chromosome 16, reference assembly (C57BL/6J)
>gi|149313536|ref|NC_000083.5|NC_000083 Mus musculus chromosome 17, reference assembly (C57BL/6J)
>gi|149321426|ref|NC_000084.5|NC_000084 Mus musculus chromosome 18, reference assembly (C57BL/6J)
>gi|149323268|ref|NC_000085.5|NC_000085 Mus musculus chromosome 19, reference assembly (C57BL/6J)
>gi|149361525|ref|NC_000086.6|NC_000086 Mus musculus chromosome X, reference assembly (C57BL/6J)
>gi|149361526|ref|NC_000087.6|NC_000087 Mus musculus chromosome Y, reference assembly (C57BL/6J) 
*/

 
 unsigned int g_unGINums[24] = {149288852,149338249,149352351,149354223,149354224,149361431,
                                149361432,149361523,149361524,149288869,149288871,149292731,
                                149292733,149292735,149301884,149304713,149313536,149321426,
                                149323268,0,0,0,149361525,149361526};

 extern char * g_geneIds[20];

 CSXReadSNV::CSXReadSNV(){
   m_pcInFile=m_pcOutFile=m_pcStatsFile=m_pcFreqFile=m_pcChroFreqFile=
   m_pcStatsFile=m_pcSNPStatsFile=m_pcINSStatsFile=m_pcDELStatsFile=
   m_pcChroStrandFile=m_pcChromOffsetFile=m_pcSNVTbl=m_pcNewFile=
   m_pcGene=m_pcGeneKeyword=m_pcExon=m_pcGeneName=m_pcdbSnp=m_pcdbIndel=
   m_pcSampleID=0;

   m_punGINums=0;

   m_bOut2Screen=m_bFQryPos=m_bFQltyScore=m_bMemMapped=m_bStateRpt=
   m_bSNVLst=m_bFilter=/*m_bSNPDensity=m_bINDELDensity=*/
   m_bDensity=m_bSNVTrace=m_bAvgQScoreLst=false; 

   m_nFQryPos1=m_nFQryPos2=m_nFQltyScore1=m_nFQltyScore2=-999;
   m_unMinSupportingReads=0; m_unMaxSupportingReads=70; m_unMaxReadDensity = 100;
   //m_fAvgQScore=
   m_fMinReadStrength=0;
   m_nQScoreOfs=33; //Default to 33
   m_eSNV = eALL;
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


inline void CSXReadSNV::PrintRemarksParam(FILE *pfParam)
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

/*
    FILE *pf = fopen("Param.ini","r+");
    char acbuf[1024]; size_t nSize=sizeof(acbuf);
    char *pChr; long lOffset=0; bool bUpdated = false;
    char actgt[1024];
    if (m_bFilter) {
        sprintf(actgt,"FRemarks=%d<=QryPos<=%d, %d<=QScore<=%d\n",
                nQryPos1,nQryPos2,nQltyScore1,nQltyScore2);
    }
    else {
        sprintf(actgt,"FRemarks=QryPos<=%d,QryPos>=%d, QScore<=%d,QScore>=%d\n",
                nQryPos1,nQryPos2,nQltyScore1,nQltyScore2);
    }

    while(!feof(pf))
    {
        if (!fgets(acbuf,nSize,pf)) break;
        if (acbuf[0] == '#') continue;

        lOffset =0-strlen(acbuf);

        pChr = strtok(acbuf,"=");

        while (pChr)
        {
            if (strcmp(pChr,"FRemarks")==0)
            {               
                fseek(pf,lOffset,SEEK_CUR);
                fwrite(actgt,strlen(actgt),1,pf);
                bUpdated = true;
                break;
            }
            pChr = strtok(NULL,"=");
        }
    }

    if (!bUpdated)
        fprintf(pf,"%s",actgt);

    if (pf) fclose(pf);
*/
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
                                   const char *pcDensFilePath,
                                   unsigned int unMinSupportingReads,
                                   unsigned int unMaxSupportingReads,
                                   float fMinReadStrength,int nMinQScore,
                                   const char *pcSampleID)
{
    m_pcdbSnp=pcdbSnp; m_pcdbIndel=pcdbIndel; m_pcDensFilePath=pcDensFilePath;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;
    m_fMinReadStrength = fMinReadStrength; m_nMinQScore=nMinQScore;
    m_pcSampleID=pcSampleID;
    m_bAvgQScoreLst=true;
}

void CSXReadSNV::setSNVLstInputs(const char *pcGene, const char *pcGeneKeyword,
                                 const char *pcExon, const char *pcCNV,
                                 unsigned int unMinSupportingReads,
                                 unsigned int unMaxSupportingReads,
                                 char cSNVType, const char *pcSampleID)
{
    m_pcGene=pcGene;m_pcGeneKeyword=pcGeneKeyword;
    m_pcExon=pcExon;m_pcCNV=pcCNV;
    m_unMinSupportingReads=unMinSupportingReads;
    m_unMaxSupportingReads=unMaxSupportingReads;    
    m_eSNV = (eSNVType)cSNVType;
    m_pcSampleID=pcSampleID;
    m_bSNVLst = true;
}

/*
void CSXReadSNV::setSNPDensityInput(const char* pcFile, const char* pcSampleID)
{
    m_pcInFile=pcFile; m_pcSampleID=pcSampleID; m_bSNPDensity = true;
}


void CSXReadSNV::setINDELDensityInput(const char* pcFile, const char* pcSampleID)
{
    m_pcInFile=pcFile; m_pcSampleID=pcSampleID; m_bINDELDensity = true;
}
*/

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
    if (m_pSNVClustered) delete[] m_pSNVClustered;
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
       if (rc == -1) {fprintf(stderr, "munmap failed"); return;}
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
    fprintf(pfStats,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfStats,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfStats,"#REMARKS\t");PrintRemarksParam(pfStats);
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
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pf,"#REMARKS\t");PrintRemarksParam(pf);
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


bool CSXReadSNV::LoadGeneKeyword(FILE *pf)
{
    if (!pf) return false;

    GENE_KEYWORD_LIST::iterator itr;
    char acbuf[516];  size_t nBufSize=sizeof(acbuf);
    char *pChr, *pcKey,*pcKeyword; unsigned int unLen;
    while (!feof(pf))
    {
        if (!fgets(acbuf,nBufSize,pf)) break;

        pChr = strchr(acbuf,'\t');
        if (!pChr) continue;
        unLen = pChr-acbuf;
        pcKey = new char[unLen+1];
        strncpy(pcKey,acbuf,unLen);
        pcKey[unLen]=0;

        unLen = strlen(acbuf) - unLen-3;//omit \t and \n
        pcKeyword = new char[unLen+1];
        strncpy(pcKeyword,pChr+1,unLen);
        pcKeyword[unLen]=0;

        itr = m_GeneKeywordList.find(pcKey);

        if (itr==m_GeneKeywordList.end())
        {
            m_GeneKeywordList[pcKey]=pcKeyword;
        }
    }
    
    return true;
}


void CSXReadSNV::OutputSNVList()
{    
    FILE *pfSrc=NULL,*pfGene=NULL,*pfGeneKeyword=NULL,*pfExon=NULL,*pfCNV=NULL;
    char *pChr=NULL;

    pfSrc = fopen(m_pcInFile,"r"); if (!pfSrc) {fprintf(stdout,"Failed to open %s ...\n",m_pcInFile); goto Exit;}
    pfGene = fopen(m_pcGene, "r"); if (!pfGene) {fprintf(stdout,"Failed to open %s ...\n",m_pcGene); goto Exit;}
    pfGeneKeyword = fopen(m_pcGeneKeyword, "r"); if (!pfGeneKeyword) {fprintf(stdout,"Failed to open %s ...\n",m_pcGeneKeyword); goto Exit;}
    pfExon = fopen(m_pcExon, "r"); if (!pfExon) {fprintf(stdout,"Failed to open %s ...\n",m_pcExon); goto Exit;}    
    pfCNV = fopen(m_pcCNV,"r");if (!pfCNV) {fprintf(stdout,"Failed to open %s ...\n",m_pcCNV); goto Exit;}
    
    fprintf(stdout,"Generating Gene Map...\n"); if (!GenerateGeneMap(pfGene)) goto Exit;
    fprintf(stdout,"Generating GeneKeyword Map...\n"); if (!LoadGeneKeyword(pfGeneKeyword)) goto Exit;
    fprintf(stdout,"Generating Exon Map...\n"); if (!GenerateExonMap(pfExon)) goto Exit;
    fprintf(stdout,"Generating CNV Map...\n"); if (!GenerateCNVMap(pfCNV)) goto Exit;
   
    pChr = strchr(m_pcGene,'7');
    m_punGINums = g_unGINums;//(!pChr)?g_unGINums_G36:g_unGINums_G37;

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

    if (pfSrc) fclose(pfSrc); if (pfGene) fclose(pfGene);
    if (pfGeneKeyword) fclose(pfGeneKeyword);
    if (pfExon) fclose(pfExon); if (pfCNV) fclose(pfCNV);

    ClrGeneMap(); ClrExonMap(); ClrCNVMap(); m_GeneKeywordList.clear();
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


void CSXReadSNV::PrintGeneKeyword(FILE *pf)
{
    int i=0; bool bFound=false;

    for (;i<20;i++)
    {
        if (!g_geneIds[i]) break;

        m_GeneKeywordItr = m_GeneKeywordList.find(g_geneIds[i]);

        if (m_GeneKeywordItr!=m_GeneKeywordList.end()){
            fprintf(pf,"%s;",m_GeneKeywordItr->second);
            if (!bFound) bFound=true; 
        }
    }

    if (!bFound) fprintf(pf,"-");    
}


void CSXReadSNV::PrintSNPList(FILE* pfSrc)
{
    char acbuf[1024],acAllele[4],acdbSnp[1024], acAvgQScore[25],acHomoHet[10];
    char *pChr=NULL; float fLocalCNV=0.0;
    float fAvgReadDens=0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulSNPReads=0,ulReadDensity=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nSNPIdx=0/*,nReadStrength=0*/; stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    
    FILE *pfSNP,*pfSyno,*pfSNPStats; pfSNP=pfSyno=pfSNPStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)];

    int nTotalTabs=0;

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
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfSNP,"#REMARKS\t");PrintRemarksParam(pfSNP);fprintf(pfSNP,"\n");
    fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"Chromosome\tGiNumber\tOffset\tNucleotide_Variant\tSNP_Reads"
                  "\tTotal_Read_Density\tAvg_Read_Density\tAvg_QScore\tdbSNP\t"
                  "Gene_Name\tGene_Description\tKeyword\tExon"
                  "\tZygosity\tLocal Copy Number"
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

        fprintf(pfSNP,"%s",pChr);//[0]);

        if (!isalpha(pChr[0])) unChromosome = atoi(pChr); 
        else if (pChr[0]=='X') unChromosome = 23;
        else if(pChr[0]=='Y') unChromosome = 24;       
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        if (nTotalTabs == 15){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulSNPReads=atol(pChr); //Total_SNP_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); strcpy(acAvgQScore,pChr);//Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); // Read Strength
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local CNV
        }
        else if(nTotalTabs == 8){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); ulSNPReads=atol(pChr); //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t"); strcpy(acAvgQScore,pChr);//Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            //pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }  

        unGINum = m_punGINums[unChromosome-1];

        /*
        fprintf(pfSNP,"\t%u\t%u\t%s\t%u\t%u\t%s\t",
                unGINum,ulOffset,acAllele,ulSNPReads,ulReadDensity,acdbSnp);
        */

        fprintf(pfSNP,"\t%u\t%u\t%s\t%u\t%u\t",
                unGINum,ulOffset,acAllele,ulSNPReads,ulReadDensity);

        fAvgReadDens = Division(ulReadDensity,ulSNPReads);
        fprintf(pfSNP,"%.2f\t%.2f\t%s\t",fAvgReadDens,atof(acAvgQScore),acdbSnp);         
 
        bGene=OutputGeneIDnDesc(unChromosome,ulOffset,pfSNP); fprintf(pfSNP,"\t");
        PrintGeneKeyword(pfSNP); fprintf(pfSNP,"\t");
        

        if (bGene){
            bExon=OutputExonID(unChromosome,ulOffset,pfSNP);

            fprintf(pfSNP,"\t");
        }
        else{
            bExon = false; fprintf(pfSNP,"-\t");
        }
        
        /*
        fAvgReadDens = Division(ulReadDensity,ulSNPReads);
        fprintf(pfSNP,"%s\t%.2f\t",acAvgQScore,fAvgReadDens);
        */

        if (nTotalTabs != 8) {
           fprintf(pfSNP,"%s\t",acHomoHet); fprintf(pfSNP,"%.2f\t",fLocalCNV);
        } 
        else fprintf(pfSNP,"-\t-\t");//,acHomoHet);  

        OutputCNVSearch(unChromosome,ulOffset,pfSNP);fprintf(pfSNP,"\t");                
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

            for (int i=0;i<20;i++) {
                if (g_geneIds[i]){
                    if (i==0) fprintf(pfSyno,"%s",g_geneIds[i]);
                    else fprintf(pfSyno,";%s",g_geneIds[i]);
                }
            }
            fprintf(pfSyno,"\n");
       }                                
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
    char acbuf[1024],acAllele[4],acdbSnp[1024], acAvgQScore[25],acHomoHet[10]; char *pChr=NULL;
    float fAvgReadDens=0.0,fLocalCNV=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulINSReads=0,ulReadDensity=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nINSIdx/*,nReadStrength=0*/; stSNVTblStats astINSStats[m_unMaxSupportingReads+1];

    FILE *pfINS,*pfINSStats; pfINS=pfINSStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)]; int nTotalTabs=0;

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
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfINS,"#REMARKS\t");PrintRemarksParam(pfINS);fprintf(pfINS,"\n");
    fprintf(pfINS,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfINS,"Chromosome\tGiNumber\tOffset\tInserted_Base\tINS_Reads"
                  "\tTotal_Read_Density\tAvg_Read_Density\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tExon"
                  "\tZygosity\tLocal Copy Number\tKnown_CNV_Region\n");

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

        fprintf(pfINS,"%s",pChr);//[0]);

        if (!isalpha(pChr[0])) unChromosome = atoi(pChr); 
        else if (pChr[0]=='X') unChromosome = 23;
        else if (pChr[0]=='Y') unChromosome = 24;        
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        if (nTotalTabs == 15){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulINSReads=atol(pChr); //Total_SNP_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); strcpy(acAvgQScore,pChr);//Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); //Read Strength
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local CNV
        }
        else if(nTotalTabs == 8){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Inserted base
            pChr = strtok(NULL,"\t"); ulINSReads=atol(pChr); //SNP_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t");  strcpy(acAvgQScore,pChr);//Avg_QScore
            pChr = strtok(NULL,"\t");  strcpy(acdbSnp,pChr); //dbSnp
            //pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];
        fprintf(pfINS,"\t%u\t%u\t%s\t%u\t%u\t",unGINum,ulOffset,acAllele,ulINSReads,ulReadDensity);

        fAvgReadDens = Division(ulReadDensity,ulINSReads);
        fprintf(pfINS,"%.2f\t%.2f\t%s\t",fAvgReadDens,atof(acAvgQScore),acdbSnp);

        bGene=OutputGeneIDnDesc(unChromosome,ulOffset,pfINS); fprintf(pfINS,"\t");
        PrintGeneKeyword(pfINS); fprintf(pfINS,"\t");

        if (bGene){
            bExon=OutputExonID(unChromosome,ulOffset,pfINS);  fprintf(pfINS,"\t");
        }
        else{
            bExon = false; fprintf(pfINS,"-\t");
        }

        if (nTotalTabs != 8) {
            fprintf(pfINS,"%s\t",acHomoHet); fprintf(pfINS,"%.2f\t",fLocalCNV);
        }
        else fprintf(pfINS,"-\t-\t");

        OutputCNVSearch(unChromosome,ulOffset,pfINS);fprintf(pfINS,"\n");

        nINSIdx = (ulINSReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulINSReads-1;

        astINSStats[nINSIdx].unCnt++;
        astINSStats[nINSIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astINSStats[nINSIdx].undbSNP++;
        if (bGene) astINSStats[nINSIdx].unGene++;
        if (bExon) astINSStats[nINSIdx].unExon++;
    }

    fprintf(stdout,"Printing INS Statistic Report...\n");

    if (pfINSStats) PrintStatsRpt(pfINSStats,astINSStats,"INS");

ExitFunc:
    if (pfINS) fclose(pfINS); if (pfINSStats) fclose(pfINSStats);
}


void CSXReadSNV::PrintDELList(FILE* pfSrc)
{
    char acbuf[1024],acAllele[4],acdbSnp[1024], acAvgQScore[25],acHomoHet[10]; char *pChr=NULL;
    float fAvgReadDens=0.0, fLocalCNV=0.0; unsigned int unChromosome=0,unGINum=0;
    unsigned long ulOffset=0,ulDELReads=0,ulReadDensity=0,ulTest1=0,ulTest2=0;
    bool bGene,bExon; int nDELIdx/*,nReadStrength=0*/; stSNVTblStats astDELStats[m_unMaxSupportingReads+1];

    FILE *pfDEL,*pfDELStats; pfDEL=pfDELStats=0;

    struct tm *ptm; time_t ttime; (void) time(&ttime); ptm = localtime(&ttime);
    char acFileName[100+strlen(m_pcSampleID)]; int nTotalTabs=0;

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
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfDEL,"#REMARKS\t");PrintRemarksParam(pfDEL);fprintf(pfDEL,"\n");
    fprintf(pfDEL,"#FILTER\t%u<=Supporting Reads<=%u, Read Strength>=%.2f\n\n\n",
                  m_unMinSupportingReads,m_unMaxSupportingReads, m_fMinReadStrength);

    fprintf(pfDEL,"Chromosome\tGiNumber\tOffset\tDeleted_Base\tTotal_DEL_Reads"
                  "\tTotal_Read_Density\tAvg_Read_Density\tAvg_QScore\tdbSNP"
                  "\tGene_Name\tGene_Description\tKeyword\tExon"
                  "\tLocal Copy Number\tZygosity\ttKnown_CNV_Region\n");

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

        fprintf(pfDEL,"%s",pChr);//[0]);

        if (!isalpha(pChr[0])) unChromosome = atoi(pChr);
        else if (pChr[0]=='X') unChromosome = 23;
        else if (pChr[0]=='Y') unChromosome = 24;        
        else {fprintf(stdout,"Invalid Chromosome number %d ...\n",atoi(pChr)); continue;}

        if (nTotalTabs == 15){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Nucleotide_Variant
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulDELReads=atol(pChr); //Total_DEL_Reads
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Total_Read_Density
            pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
            pChr = strtok(NULL,"\t"); strcpy(acAvgQScore,pChr);//Total_Avg_QScore
            pChr = strtok(NULL,"\t"); strcpy(acdbSnp,pChr); //dbSnp
            pChr = strtok(NULL,"\t"); //Read Strength
            pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
            pChr = strtok(NULL,"\t"); fLocalCNV = atof(pChr); //Local CNV
        }
        else if(nTotalTabs == 7){
            pChr = strtok(NULL,"\t"); ulOffset = atol(pChr); //Offset
            pChr = strtok(NULL,"\t"); strcpy(acAllele,pChr); //Deleted_Variant
            pChr = strtok(NULL,"\t"); ulDELReads=atol(pChr); //DEL_Reads
            pChr = strtok(NULL,"\t"); ulReadDensity=atol(pChr); //Read_Density
            pChr = strtok(NULL,"\t");  strcpy(acAvgQScore,pChr);//Avg_QScore
            pChr = strtok(NULL,"\t");  strcpy(acdbSnp,pChr); //dbSnp
            //pChr = strtok(NULL,"\t"); strcpy(acHomoHet,pChr); //Homo_Het
        }
        else {
           fprintf(stdout,"Number of columns in input file is not correct...\n");
           fprintf(stdout,"It has %d tabs\n",nTotalTabs);
           return;
        }

        unGINum = m_punGINums[unChromosome-1];

        fprintf(pfDEL,"\t%u\t%u\t%s\t%u\t%u\t",unGINum,ulOffset,acAllele,ulDELReads,ulReadDensity);

        fAvgReadDens = Division(ulReadDensity,ulDELReads);
        fprintf(pfDEL,"%.2f\t%.2f\t%s\t",fAvgReadDens,atof(acAvgQScore),acdbSnp);  
         
        bGene=OutputGeneIDnDesc(unChromosome,ulOffset,pfDEL); fprintf(pfDEL,"\t");
        PrintGeneKeyword(pfDEL); fprintf(pfDEL,"\t");

        if (bGene){
            bExon=OutputExonID(unChromosome,ulOffset,pfDEL);

            fprintf(pfDEL,"\t");
        }
        else{
            bExon = false; fprintf(pfDEL,"-\t");
        }

        if (nTotalTabs != 7) {
            fprintf(pfDEL,"%s\t",acHomoHet); fprintf(pfDEL,"%.2f\t",fLocalCNV);
        }
        else fprintf(pfDEL,"-\t-\t");

        OutputCNVSearch(unChromosome,ulOffset,pfDEL);fprintf(pfDEL,"\n");

        nDELIdx = (ulDELReads > m_unMaxSupportingReads)?m_unMaxSupportingReads:ulDELReads-1;

        astDELStats[nDELIdx].unCnt++;
        astDELStats[nDELIdx].ulTotalReadDens += ulReadDensity;

        if (acdbSnp[0]!='-') astDELStats[nDELIdx].undbSNP++;
        if (bGene) astDELStats[nDELIdx].unGene++;
        if (bExon) astDELStats[nDELIdx].unExon++;
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
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pf,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");    
    fprintf(pf,"#REMARKS\t");PrintRemarksParam(pf);fprintf(pf,"\n");
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
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pf,"#REMARKS\t");PrintRemarksParam(pf);
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


inline unsigned int CSXReadSNV::GetReadDensityEx(FILE *pf, unsigned int uioffset, unsigned int unFranking)
{
    if (!pf) return 0;

    unsigned int unRead=0, unTotalReads=0, unFinalIdx=unFranking+1;
    uioffset--;

    for (unsigned int i=unFranking;i<unFinalIdx;i++)//for (int i=1;i<2;i++)
    {
         fseek(pf,(uioffset-i)*2,SEEK_SET);
         fread(g_acBuf,2,1,pf);

         unRead = g_acBuf[1];unRead<<=8;
         unRead|= g_acBuf[0];
         unTotalReads += unRead;
    }

    for (unsigned int i=unFranking;i<unFinalIdx;i++)//for (int i=1;i<2;i++)
    {
         fseek(pf,(uioffset+i)*2,SEEK_SET);
         fread(g_acBuf,2,1,pf);

         unRead = g_acBuf[1];unRead<<=8;
         unRead|= g_acBuf[0];
         unTotalReads += unRead;
    }
 
    g_unReads = (unsigned int)((double)((double)(unTotalReads)/(unFranking*2))+0.5);

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
    if (((const stSNVRec *) a)->cChromosome != ((const stSNVRec *) b)->cChromosome)
       return ((const stSNVRec *) a)->cChromosome - ((const stSNVRec *) b)->cChromosome;

    if (((const stSNVRec *) a)->unOffset != ((const stSNVRec *) b)->unOffset)
        return ((const stSNVRec *) a)->unOffset - ((const stSNVRec *) b)->unOffset;

    if (((const stSNVRec *) a)->cRefBase != ((const stSNVRec *) b)->cRefBase)
        return ((const stSNVRec *) a)->cRefBase - ((const stSNVRec *) b)->cRefBase;

    return ((const stSNVRec *) a)->cVarBase - ((const stSNVRec *) b)->cVarBase;
}


void CSXReadSNV::AllocMergeListSize()
{
    //fprintf(stdout,"------------------------------------------------\n");
    //fprintf(stdout,"INIT_MERGED_LIST_SIZE = %u\n",INIT_MERGED_LIST_SIZE);
    //fprintf(stdout,"EXT_MERGED_LIST_SIZE = %u\n",INIT_MERGED_LIST_SIZE);
    //fprintf(stdout,"Initial merge_list_alloc = %u\n",merge_list_alloc);
    //fprintf(stdout,"merge_list_used = %u\n",merge_list_used);

	if (m_pSNVClustered == NULL) {
		merge_list_alloc = INIT_MERGED_LIST_SIZE;
		merge_list_used = 0;
		try {
			m_pSNVClustered = new stSNVClusterRec[merge_list_alloc];
                        //fprintf(stdout,"Initial Allocation\n");
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
		stSNVClusterRec *pSNVClustered_tmp = new stSNVClusterRec[merge_list_alloc];
                memcpy(pSNVClustered_tmp,m_pSNVClustered,sizeof(stSNVClusterRec)*(merge_list_alloc - EXT_MERGED_LIST_SIZE));
                delete [] m_pSNVClustered;
		m_pSNVClustered = pSNVClustered_tmp;
                //fprintf(stdout,"Addition Allocation - merge_list_alloc = %u\n",merge_list_alloc);
	}
	catch (...) {
		fprintf(stderr, "Error: insufficient memory for Pair Info list\n");
		fprintf(stderr, "Program terminated.\n");
		exit(1);
	}

        //fprintf(stdout,"------------------------------------------------\n");
}


inline void CSXReadSNV::GenerateSNVTbl()
{
     
    FILE *pfIn=fopen("unsort.txt","w");

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
           fprintf(pfIn,"%d\t%u\t%c\t%c\t%d\t%d\t%d\n",
                m_pSNV[i].cChromosome,ComputeOffset(m_pSNV[i].acOffset),
                toupper(m_pSNV[i].cRefBase),toupper(m_pSNV[i].cVarBase),
                m_pSNV[i].ucQltyScore,m_pSNV[i].ucQryPos,m_pSNV[i].sPEnd);
    }
    fclose(pfIn);
   

    unsigned long ulTest1=0, ulTest2=0;

#if _TESTING
    m_ulTotalRecs=25000; //Testing
#endif

    stSNVRec *pTmp = new stSNVRec[m_ulTotalRecs];

    fprintf(stdout,"Converting chromosome to absolute value and compute offset...\n");

    char cRefBase, cVarBase; char cChromosome;

    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        cRefBase=toupper(m_pSNV[i].cRefBase);
        cVarBase=toupper(m_pSNV[i].cVarBase);

        cChromosome = abs(m_pSNV[i].cChromosome);
        pTmp[i].cChromosome = (cChromosome < 86)?cChromosome-66:cChromosome-63;
        pTmp[i].unOffset = ComputeOffset(m_pSNV[i].acOffset);
        pTmp[i].cRefBase = cRefBase;
        pTmp[i].cVarBase = cVarBase;
        pTmp[i].ucQryPos = m_pSNV[i].ucQryPos;
        pTmp[i].ucQltyScore = m_pSNV[i].ucQltyScore;
        pTmp[i].sPEnd = m_pSNV[i].sPEnd;
    }

    ClrMemoryMap();

    fprintf(stdout,"sorting...\n");
    qsort(pTmp,m_ulTotalRecs,sizeof(stSNVRec), CompareSNVRecs);

         
    FILE *pf=fopen("sorted.txt","w");
    for (unsigned long i=0; i<m_ulTotalRecs; i++)
    {
        fprintf(pf,
                "%d\t%u\t%c\t%c\t%d\t%d\t%d\n",pTmp[i].cChromosome,pTmp[i].unOffset,
                toupper(pTmp[i].cRefBase),toupper(pTmp[i].cVarBase),
                pTmp[i].ucQltyScore,pTmp[i].ucQryPos,pTmp[i].sPEnd);

        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }

#if _TESTING

        if (i+1 == 25000)
        {
            fprintf(stdout,"Total Recs Processed = %u\n\n\n",i+1);
            break;
         }
#endif
    }

    fclose(pf);  
      
    fprintf(stdout,"alloc memory for clustering...\n");

    AllocMergeListSize();

    fprintf(stdout,"clustering...\n");

    //FILE *pfLog = fopen("InvalidRec.lst","w");
    char cChromosome_Curr,cRefBase_Curr,cVarBase_Curr;
    char cRef_ACGT,cVar_ACGT;
    unsigned int unOffset_Curr,unOffset;
    unsigned short usQS; short sPEnd;

    sPEnd = pTmp[0].sPEnd;
    usQS = pTmp[0].ucQltyScore-m_nQScoreOfs;
    cChromosome_Curr = pTmp[0].cChromosome;
    unOffset_Curr = pTmp[0].unOffset;
    cRefBase_Curr = pTmp[0].cRefBase;
    cVarBase_Curr = pTmp[0].cVarBase;

    m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
    m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
    m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
    m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
    m_pSNVClustered[merge_list_used].ulTotalQS = usQS;
    m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
    m_pSNVClustered[merge_list_used].sMaxPEnd = sPEnd;
    m_pSNVClustered[merge_list_used].unCnt = 1;

    if ((pTmp[0].cRefBase!='A' && pTmp[0].cRefBase!='C' && pTmp[0].cRefBase!='G' &&
         pTmp[0].cRefBase!='T' && pTmp[0].cRefBase!='-')||
        (pTmp[0].cVarBase!='A' && pTmp[0].cVarBase!='C' && pTmp[0].cVarBase!='G' &&
         pTmp[0].cVarBase!='T' && pTmp[0].cVarBase!='T'))
    {
        m_pSNVClustered[merge_list_used].bKeep=false;
        //LogInvalidRec(pfLog,pTmp[0]);
    }

    for (unsigned long i=1; i<m_ulTotalRecs; i++)
    {
        if ((pTmp[i].cRefBase!='A' && pTmp[i].cRefBase!='C' && pTmp[i].cRefBase!='G'&&
             pTmp[i].cRefBase!='T' && pTmp[i].cRefBase!='-')||
            (pTmp[i].cVarBase!='A' && pTmp[i].cVarBase!='C' &&
             pTmp[i].cVarBase!='G' && pTmp[i].cVarBase!='T' && pTmp[i].cVarBase!='-'))
        {
           // LogInvalidRec(pfLog,pTmp[i]); 
           continue;
        }

        cRef_ACGT = pTmp[i].cRefBase;
        cVar_ACGT = pTmp[i].cVarBase;
        unOffset = pTmp[i].unOffset;
        usQS = pTmp[i].ucQltyScore-m_nQScoreOfs;
        sPEnd = pTmp[i].sPEnd;

        if (cChromosome_Curr != pTmp[i].cChromosome||unOffset_Curr != unOffset||
            cRefBase_Curr != cRef_ACGT||cVarBase_Curr != cVar_ACGT)
        {
            cChromosome_Curr = pTmp[i].cChromosome;
            unOffset_Curr = unOffset;
            cRefBase_Curr = cRef_ACGT;
            cVarBase_Curr = cVar_ACGT;

            merge_list_used++;

            if (merge_list_used == merge_list_alloc)
                AllocMergeListSize();

            m_pSNVClustered[merge_list_used].cChromosome = cChromosome_Curr;
            m_pSNVClustered[merge_list_used].unOffset = unOffset_Curr;
            m_pSNVClustered[merge_list_used].cRefBase = cRefBase_Curr;
            m_pSNVClustered[merge_list_used].cVarBase = cVarBase_Curr;
            m_pSNVClustered[merge_list_used].ulTotalQS = usQS;
            m_pSNVClustered[merge_list_used].ucMaxQS = usQS;
            m_pSNVClustered[merge_list_used].sMaxPEnd = sPEnd;
            m_pSNVClustered[merge_list_used].unCnt = 1;
        }
        else
        {
            m_pSNVClustered[merge_list_used].ulTotalQS += usQS;

            if (m_pSNVClustered[merge_list_used].ucMaxQS < usQS)
                m_pSNVClustered[merge_list_used].ucMaxQS = usQS;

            if (m_pSNVClustered[merge_list_used].sMaxPEnd < sPEnd)
                m_pSNVClustered[merge_list_used].sMaxPEnd = sPEnd;

            m_pSNVClustered[merge_list_used].unCnt++;
        }

        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs gone through clustering = %uM\n", ulTest2);
            ulTest1=0;
        }
    }

    delete[] pTmp;

    //if (pfLog) fclose(pfLog);
    fprintf(stdout,"Total Clustered Recs = %u\n", merge_list_used);
    
     
    pf=fopen("clustered.txt","w");

    for (unsigned long i=0; i < merge_list_used; i++)
    {
        fprintf(pf,"%d\t%u\t%c>%c\t%u\t%u\t%u\n",
                m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                m_pSNVClustered[i].cRefBase,m_pSNVClustered[i].cVarBase,
                m_pSNVClustered[i].unCnt,m_pSNVClustered[i].ulTotalQS,
                m_pSNVClustered[i].ucMaxQS);
    }

    fclose(pf);
     
}

void CSXReadSNV::LogInvalidRec(FILE *pf,stSNVRec &Rec)
{
    fprintf(pf,"%d\t%u\t%c\t%c\t%d\t%d\n",Rec.cChromosome, Rec.unOffset,
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


void CSXReadSNV::OutputAvgQScoreTable()
{
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
        if (i > 18 && i < 22) continue;

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
    fprintf(stdout,"Filtering SNV Table...\n"); FilterSNVTbl();
    fprintf(stdout,"Printing SNP, INS, DEL Lists...\n"); PrintAvgQScoreList();
   
Exit:

    if (pfdbSnp) fclose(pfdbSnp); if (pfdbIndel) fclose(pfdbIndel);
    ClrDBSNPMap(); ClrDBINDELMap();

    for (int j=0;j<24; j++){
        if (m_apfRD[j]) fclose(m_apfRD[j]);
    }

    delete[] m_apfRD;
}




void CSXReadSNV::PrintAvgQScoreList()
{

/////////////////////////////

    stSNVTblStats astSNPStats[m_unMaxSupportingReads+1];
    stSNVTblStats astINSStats[m_unMaxSupportingReads+1];
    stSNVTblStats astDELStats[m_unMaxSupportingReads+1];
    int nSNPIdx,nINSIdx,nDELIdx;

    unsigned long ulTest1=0, ulTest2=0;
    unsigned int unChroIdx,unReadDensity=0;

    bool bdbSNP; char acRsid[20];

    float fAvgQScore;//,fAvgReadDens;

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
    fprintf(pfSNP,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfSNP,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfSNP,"#REMARKS\t");PrintRemarksParam(pfSNP);
    fprintf(pfSNP,"#FILTER\t>=1\n\n\n");
    //fprintf(pfSNP,"#FILTER\t%u<=Supporting Reads<=%u\n\n\n",
    //              m_unMinSupportingReads,m_unMaxSupportingReads);

    fprintf(pfSNP,"Chromosome\tOffset\tNucleotide_Variant\tSNP_Reads"
                  "\tRead_Density\tAvg_QScore\tdbSNP\n");

    fprintf(pfINS,"#REPORT NAME\tINS Avg. Score Distribution\n");
    fprintf(pfINS,"#PROJECT NAME\n");
    fprintf(pfINS,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfINS,"#LANE NO\tALL\n");
    fprintf(pfINS,"#GENERATED AT\t"); PrintRptDateTime(pfINS);
    fprintf(pfINS,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfINS,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfINS,"#REMARKS\t");PrintRemarksParam(pfINS);
    fprintf(pfINS,"#FILTER\t>=1\n\n\n");

    fprintf(pfINS,"Chromosome\tOffset\tInserted_Base\tINS_Reads"
                  "\tRead_Density\tAvg_QScore\tdbSNP\n");

    fprintf(pfDEL,"#REPORT NAME\tDEL Avg. Score Distribution\n");
    fprintf(pfDEL,"#PROJECT NAME\n");
    fprintf(pfDEL,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pfDEL,"#LANE NO\tALL\n");
    fprintf(pfDEL,"#GENERATED AT\t"); PrintRptDateTime(pfDEL);
    fprintf(pfDEL,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
    fprintf(pfDEL,"#REFERENCE DATABASE\tNCBI Human Genome V37.1\n");
    fprintf(pfDEL,"#REMARKS\t");PrintRemarksParam(pfDEL);
    fprintf(pfDEL,"#FILTER\t>=1\n\n\n");

    fprintf(pfDEL,"Chromosome\tOffset\tDeleted_Base\tDEL_Reads"
                  "\tRead_Density\tAvg_QScore\tdbSNP\n");

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

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfSNP,m_pSNVClustered[i].cChromosome);
            fprintf(pfSNP,"\t%u\t%c>%c\t%u",m_pSNVClustered[i].unOffset,
                    m_pSNVClustered[i].cRefBase,m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unCnt);

            fprintf(pfSNP,"\t%u\t%.2f\t%s\t",unReadDensity,fAvgQScore,acRsid);
            fprintf(pfSNP,"\n");

            nSNPIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astSNPStats[nSNPIdx].unCnt++;
            if (bdbSNP) astSNPStats[nSNPIdx].undbSNP++;
            astSNPStats[nSNPIdx].ulTotalReadDens += unReadDensity;

        }
        else if (m_pSNVClustered[i].cRefBase =='-' )
        {
            bdbSNP=OutputINDELNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                    m_pSNVClustered[i].cVarBase, pfINS,acRsid);

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;            

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensityEx(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset,1);
             
            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfINS,m_pSNVClustered[i].cChromosome);

            fprintf(pfINS,"\t%u\t%c\t%u",m_pSNVClustered[i].unOffset,
                    m_pSNVClustered[i].cVarBase,m_pSNVClustered[i].unCnt);

            fprintf(pfINS,"\t%u\t%.2f\t%s\t",unReadDensity,fAvgQScore,acRsid);
            
            fprintf(pfINS,"\n");

            nINSIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astINSStats[nINSIdx].unCnt++;
            if (bdbSNP) astINSStats[nINSIdx].undbSNP++;
            astINSStats[nINSIdx].ulTotalReadDens += unReadDensity;
        }
        else
        {
            bdbSNP=OutputINDELNovel(m_pSNVClustered[i].cChromosome,m_pSNVClustered[i].unOffset,
                                    m_pSNVClustered[i].cRefBase,pfDEL,acRsid);

            if (!m_pSNVClustered[i].bKeep && !bdbSNP) continue;

            unChroIdx=m_pSNVClustered[i].cChromosome-1;
            fAvgQScore = CalcAvgQS(m_pSNVClustered[i].ulTotalQS,m_pSNVClustered[i].unCnt);
            unReadDensity = GetReadDensity(m_apfRD[unChroIdx], m_pSNVClustered[i].unOffset)+m_pSNVClustered[i].unCnt;

            if (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads && !bdbSNP) continue;
            if (Percentage(m_pSNVClustered[i].unCnt,unReadDensity) < m_fMinReadStrength && !bdbSNP) continue;

            OutputChrom(pfDEL,m_pSNVClustered[i].cChromosome);
            fprintf(pfDEL,"\t%u\t%c\t%u",m_pSNVClustered[i].unOffset,
                    m_pSNVClustered[i].cRefBase,m_pSNVClustered[i].unCnt);
           
            fprintf(pfDEL,"\t%u\t%.2f\t%s\t",unReadDensity,fAvgQScore,acRsid);

            fprintf(pfDEL,"\n");

            nDELIdx = (m_pSNVClustered[i].unCnt > m_unMaxSupportingReads)?m_unMaxSupportingReads:m_pSNVClustered[i].unCnt-1;
            astDELStats[nDELIdx].unCnt++;
            if (bdbSNP) astDELStats[nDELIdx].undbSNP++;
            astDELStats[nDELIdx].ulTotalReadDens += unReadDensity;
        }
    }


    fprintf(stderr,"Printing AvgQScore Statistic Report...\n");

    PrintAvgQSStatsRpt(pfSNPStats,astSNPStats,"SNP");
    PrintAvgQSStatsRpt(pfINSStats,astINSStats,"INS");
    PrintAvgQSStatsRpt(pfDELStats,astDELStats,"DEL");

ExitFunc:
    if (pfSNP) fclose(pfSNP); if (pfINS) fclose(pfINS);
    if (pfDEL) fclose(pfDEL); if (pfSNPStats) fclose(pfSNPStats);
    if (pfINSStats) fclose(pfINSStats); if (pfDELStats) fclose(pfDELStats);

////////////
    
}



void CSXReadSNV::PrintAvgQSStatsRpt(FILE *pf, stSNVTblStats *pstTblStats, char *pcType)
{   
    fprintf(pf,"#REPORT NAME\t%s Avg. QScore Statistic\n",pcType);
    fprintf(pf,"#PROJECT NAME\n");
    fprintf(pf,"#SAMPLE ID\t%s\n",m_pcSampleID);
    fprintf(pf,"#LANE NO\tALL\n");
    fprintf(pf,"#GENERATED AT\t"); PrintRptDateTime(pf);
    fprintf(pf,"#PROGRAM & BUILD\tReadSNV Rev 4.0.0\n");
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








