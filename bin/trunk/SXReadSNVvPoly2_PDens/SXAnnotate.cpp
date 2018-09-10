#include "SXAnnotate.h"


SXAnnotate::SXAnnotate()
{
}


SXAnnotate::~SXAnnotate()
{
    if(m_snpfile.is_open()) m_snpfile.close();
    if (m_pcAnnotate) delete [] m_pcAnnotate;       
}


bool SXAnnotate::LoadDataSource(const char* pcbinnedfile, 
                                const char* pcExceptionfilename,
                                const char* pcfilename)
{
    m_snpfile.open(pcfilename, ios::in | ios::binary);

    if(!m_snpfile.is_open()){
        cout<<"could not open " << pcfilename << "...\n"; return false;
    }

    m_nfdAnnotate =open(pcfilename, O_RDONLY, S_IRUSR);
    if (m_nfdAnnotate == -1) {
        fprintf(stdout, "File %s not found...\n",pcfilename); return false;
    }

    struct stat st; fstat(m_nfdAnnotate, &st); m_file_size = (size_t)st.st_size;

    char *pcRecs=(char*)mmap(NULL,m_file_size,PROT_READ,MAP_PRIVATE,m_nfdAnnotate,0);

    if (pcRecs == (void *)-1) {
        fprintf(stdout, "Error: mmap failed: %s\n", strerror(errno));return false;
    }

    m_pcAnnotate = new char[m_file_size+1]; strcpy(m_pcAnnotate,pcRecs);

    int rc = munmap((char *)pcRecs, m_file_size);
    if (rc == -1) fprintf(stdout, "munmap failed...\n");
    close(m_nfdAnnotate);

    if (!LoadExceptionfile(pcExceptionfilename)) return false;

    return Loadbinfile(pcbinnedfile);
}


extern "C" int CompareBinnedData(const void *a, const void *b)
{
    return ((const stBinnedData *)a)->nOffset - ((const stBinnedData *)b)->nOffset;
}

bool SXAnnotate::Loadbinfile(const char* pcfilename)
{

    FILE *pf = fopen(pcfilename,"r");

    if (!pf){fprintf(stdout, "File %s not found...\n",pcfilename); return false;}

    fprintf(stdout,"Loading Binned file...\n");

    char line[40960]; int nLen=sizeof(line);    
    char *pChr; stBinnedIdx *pBinnedIdx=NULL; 

    while(!feof(pf))
    {
        if(!fgets(line,nLen,pf)) break;
        
        pChr = strtok(line,"\t");

        if (!isalpha(pChr[0])) pBinnedIdx = &m_aBinnedIdx[atoi(pChr)-1];        
        else if (pChr[0]=='x') pBinnedIdx = &m_aBinnedIdx[22];
        else pBinnedIdx = &m_aBinnedIdx[23];
        
        if (pBinnedIdx->nUsedSubscript == pBinnedIdx->nTotalSubscript){
            pBinnedIdx->nTotalSubscript+=REC_SIZE;
            pBinnedIdx->pBinnedData = (stBinnedData*)realloc(pBinnedIdx->pBinnedData,
                                       sizeof(stBinnedData)*pBinnedIdx->nTotalSubscript);
        }

        pChr = strtok(NULL,"\t");         
        pBinnedIdx->pBinnedData[pBinnedIdx->nUsedSubscript].nOffset = atoi(pChr);
        pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t");
        pBinnedIdx->pBinnedData[pBinnedIdx->nUsedSubscript].nIdx = atoi(pChr);
        pBinnedIdx->nUsedSubscript++;
    }

    fclose(pf);

    for (int i=0;i<24;i++){
        qsort(m_aBinnedIdx[i].pBinnedData,m_aBinnedIdx[i].nUsedSubscript,
              sizeof(stBinnedData), CompareBinnedData);
    }
    return true;
}


bool SXAnnotate::LoadExceptionfile(const char* pcfilename)
{
    FILE *pf = fopen(pcfilename,"r");

    if (!pf){fprintf(stdout, "File %s not found...\n",pcfilename); return false;}

    fprintf(stdout,"Loading Exception file...\n");

    char line[40960], acData[40960],*pChr; int nLen=sizeof(line);
    stExceptionIdx *pExceptionIdx=NULL;

    while(!feof(pf))
    {
        if(!fgets(line,nLen,pf)) break;

        pChr = strchr(line,'\r'); if (pChr) *pChr=0;
        pChr = strchr(line,'\n'); if (pChr) *pChr=0;

        strcpy(acData,line); pChr = strtok(line,"\t");

        if (!isalpha(pChr[0])) pExceptionIdx = &m_aExceptionIdx[atoi(pChr)-1];
        else if (pChr[0]=='x') pExceptionIdx = &m_aExceptionIdx[22];
        else pExceptionIdx = &m_aExceptionIdx[23];

        if (pExceptionIdx->nUsedSubscript == pExceptionIdx->nTotalSubscript){
            pExceptionIdx->nTotalSubscript+=REC_SIZE;
            pExceptionIdx->pExceptionData = (stExceptionData*)realloc(pExceptionIdx->pExceptionData,
                                             sizeof(stExceptionData)*pExceptionIdx->nTotalSubscript);
        }
        pChr = strtok(NULL,"\t");
        pExceptionIdx->pExceptionData[pExceptionIdx->nUsedSubscript].nStart = atoi(pChr);
        pChr = strtok(NULL,"\t"); 
        pExceptionIdx->pExceptionData[pExceptionIdx->nUsedSubscript].nStop = atoi(pChr);
        pExceptionIdx->pExceptionData[pExceptionIdx->nUsedSubscript].pcData = new char[strlen(acData)+1];
        strcpy(pExceptionIdx->pExceptionData[pExceptionIdx->nUsedSubscript].pcData,acData);
        pExceptionIdx->nUsedSubscript++;
    }

    fclose(pf);  return true;
}

  
void SXAnnotate::SearchRecords(size_t index,char cChrnum,size_t qstart, size_t qstop)
{    
    char *pcTab,*pcRec,*pChr;  size_t recstart,recstop; int i;
   
    while(m_pcAnnotate[index])
    {
        pcRec = &m_pcAnnotate[index]; i=0;

        do{
            index++; i++;
        }while (m_pcAnnotate[index] != '\n');

        index++; i++;  
 
        if (pcRec[1]=='\t'){
            if (pcRec[0] < 65){
                if ((pcRec[0]-48)!= cChrnum) return;
            }
            else{
                if (pcRec[0]!= cChrnum) return;
            }
            pChr = &pcRec[2];
        }
        else{
            pcRec[2]='\0';
            if (atoi(pcRec)!=cChrnum) return;
            pcRec[2]='\t';
            pChr = &pcRec[3];
        }

        pcTab = strchr(pChr,'\t'); *pcTab='\n';
        recstart = atoi(pChr); *pcTab='\t';

        if (qstart < recstart) return;

        pChr = pcTab+1; pcTab = strchr(pChr,'\t'); *pcTab='\n';
        recstop = atoi(pChr); *pcTab='\t';
        
        if ((recstop-((recstart/50000)*50000)) > 50000) continue;

        if (!(recstart > qstop) && !(recstop < qstart)){ 
            if (m_Records.nTotalSize < (m_Records.nCurrSize+i+1)) // 1 for null
            {   m_Records.nTotalSize += REC_SIZE;
                m_Records.pcData = (char*)realloc(m_Records.pcData,m_Records.nTotalSize);
            }
            memcpy(m_Records.pcData+m_Records.nCurrSize,pcRec,i);
            m_Records.nCurrSize+=i;
            m_Records.pcData[m_Records.nCurrSize] = '\0';                        
        }
    }
    
    return;
} 


int SXAnnotate::GetIndex(char cChromosome, int nOffset)
{
    stBinnedIdx *pBinnedIdx = NULL;

    if (!isalpha(cChromosome))pBinnedIdx = &m_aBinnedIdx[cChromosome-1];
    else if (cChromosome=='x') pBinnedIdx = &m_aBinnedIdx[22];
    else pBinnedIdx = &m_aBinnedIdx[23];

    int nStart,nLast, nSeek;
    nStart=0; nLast=pBinnedIdx->nUsedSubscript;

    do
    {
        nSeek = (nStart+nLast)/2;

        if  (pBinnedIdx->pBinnedData[nSeek].nOffset == nOffset)
            return pBinnedIdx->pBinnedData[nSeek].nIdx;
        else if (pBinnedIdx->pBinnedData[nSeek].nOffset > nOffset)
            nLast = nSeek -1;
        else
            nStart = nSeek+1;
    }while(!(nStart > nLast));

    return -1;
}


char* SXAnnotate::GetRecord(char cChr, size_t qstart, size_t qstop)
{                     
    size_t qtemp; qtemp = qstart; 
    m_Records.pcData[0]='\0'; m_Records.nCurrSize = 0;

    int nOffset,nIdx; char cChromosome = tolower(cChr);

    while(qtemp <= qstop)
    {
        nOffset = (qtemp/50000)*50000; nIdx = GetIndex(cChromosome,nOffset);
        if (nIdx > -1) 
            SearchRecords(nIdx,cChromosome,qstart,qstop);
        
        qtemp += 50000;
    }

    GetExceptionRecord(cChr, qstart,qstop);
  
    return m_Records.pcData;
}


void SXAnnotate::GetExceptionRecord(char cChar, size_t qstart, size_t qstop)
{
    stExceptionIdx *pExceptionIdx=NULL;    

    if (cChar=='x') pExceptionIdx = &m_aExceptionIdx[22];
    else if (cChar=='y') pExceptionIdx = &m_aExceptionIdx[23];
    else pExceptionIdx = &m_aExceptionIdx[cChar-1];

    int nLen=0;
    for (int i=0; i < pExceptionIdx->nUsedSubscript; i++)
    {   
        if (qstop < pExceptionIdx->pExceptionData[i].nStart) break;

        if (!(qstart > pExceptionIdx->pExceptionData[i].nStop) &&
            !(qstop < pExceptionIdx->pExceptionData[i].nStart))
        {
            nLen = strlen(pExceptionIdx->pExceptionData[i].pcData);
            if (m_Records.nTotalSize < (m_Records.nCurrSize+nLen+1)) // 1 for null
            {
                m_Records.nTotalSize += REC_SIZE;
                m_Records.pcData = (char*)realloc(m_Records.pcData,m_Records.nTotalSize);
            }
            memcpy(m_Records.pcData+m_Records.nCurrSize,pExceptionIdx->pExceptionData[i].pcData,nLen);
            m_Records.nCurrSize+=nLen;
            m_Records.pcData[m_Records.nCurrSize] = '\n';
            m_Records.pcData[++m_Records.nCurrSize] = '\0';
        }
    }
    
    return;
}

/*
int main(int argc, char* argv[])
{
    if (argc !=5)
    {
       SXAnnotate sobj;

       sobj.LoadDataSource(argv[1],argv[2],argv[3]);

       char *pcData = sobj.GetRecord(1,51325,51325); //1942,815700);

       if (!pcData[0]) return 0;

       char *pChr=NULL;

       pChr = strtok(pcData,"\n");

       while(pChr){
           cout << pChr<<endl;pChr = strtok(NULL,"\n");
       }
        
    }
    else
    {
       FILE *pfIn = fopen(argv[4],"r");
       if (!pfIn){
           fprintf(stdout,"Failed to open %s ...\n",argv[3]); exit(9);
       }

       SXAnnotate sobj;
       
       sobj.LoadDataSource(argv[1],argv[2],argv[3]);

       char acbuf[1024],*pChr,*pChr2; char *pcData=NULL;
       int nOffset; char cChromosome;

       while(!feof(pfIn))
       {
           if (!fgets(acbuf,sizeof(acbuf),pfIn)) break;
           pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
           pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

           pChr = strtok(acbuf,"\t");

           if (!isalpha(pChr[0])) cChromosome = atoi(pChr);
           else cChromosome = pChr[0];

           pChr2 = strtok(NULL,"\t"); nOffset = atoi(pChr2);

           cout << "****" << nOffset << endl;    
                      
           pcData = sobj.GetRecord(cChromosome,nOffset,nOffset);

           if (pcData[0]){
               pChr = strtok(pcData,"\n");

               while(pChr){
                   cout << pChr<<endl;pChr = strtok(NULL,"\n");
               }
           }
       }
    }
    
  return 0;
}
*/



