#include "SXChkSeqSVar.h"


static void banner(char *argv[])
{
    printf("Synamatix SXChkSeqSVar Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: \n", argv[0]);
    printf("\t<Seq filename> <SNP list filename> <INS list filename> <DEL list filename> <Output filename>\n\n");
}


extern "C" int CompareSNVRecs(const void *a, const void *b)
{
    if (((const stRec *) a)->nOffset != ((const stRec *) b)->nOffset)
       return ((const stRec *) a)->nOffset - ((const stRec *) b)->nOffset;

    if (((const stRec *) a)->acAllele[2]!='0')
    { 
        if (((const stRec *) a)->acAllele[2] != ((const stRec *) b)->acAllele[2])
            return ((const stRec *) a)->acAllele[2] - ((const stRec *) b)->acAllele[2];
    }
    else
    {
        if (((const stRec *) a)->acAllele[0] != ((const stRec *) b)->acAllele[0])
            return ((const stRec *) a)->acAllele[0] - ((const stRec *) b)->acAllele[0]; 
    } 

}


void GenerateSNVList(int nSNVType)
{
    char acbuf[40960], *pChr=NULL,*pcLast=NULL; stRec *pstRec; 
    int nUsed, nIdx; unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while (!feof(g_stSNVList[nSNVType].pf))
    {
         if (++ulTest1 == 100000)
         {
            fprintf(stdout,"Total Recs Processed[%d] = %u x 0.1M\n", nSNVType,++ulTest2);
            ulTest1=0;
         }

         if (!fgets(acbuf,sizeof(acbuf),g_stSNVList[nSNVType].pf)) break;    
         
         pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
         pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

         if (acbuf[0]==0||acbuf[0]=='#'||acbuf[0]=='C') continue;
       
         pChr = strtok(acbuf,"\t"); 

         if (pChr[0] == 'X') nIdx = 22;
         else if (pChr[0] == 'Y') nIdx = 23;
         else nIdx = atoi(pChr)-1;

         if (g_stSNVList[nSNVType].anUsed[nIdx] == g_stSNVList[nSNVType].anTotal[nIdx])
         {
             g_stSNVList[nSNVType].anTotal[nIdx] += REC_SIZE;
             g_stSNVList[nSNVType].apstRec[nIdx] = (stRec*)realloc(g_stSNVList[nSNVType].apstRec[nIdx], 
                                                                  sizeof(stRec)*g_stSNVList[nSNVType].anTotal[nIdx]);
         }

         pstRec = g_stSNVList[nSNVType].apstRec[nIdx]; nUsed = g_stSNVList[nSNVType].anUsed[nIdx];  
         pChr = strtok(NULL,"\t"); pChr = strtok(NULL,"\t"); pstRec[nUsed].nOffset=atoi(pChr);
         pChr = strtok(NULL,"\t"); strcpy(pstRec[nUsed].acAllele,pChr);

         do {
               pcLast = pChr; pChr = strtok(NULL,"\t");
         } while(pChr);

         pstRec[nUsed].fCfdScore = atof(pcLast); g_stSNVList[nSNVType].anUsed[nIdx]++;
    }

    for (int i=0; i<24;i++)
    {
        qsort( g_stSNVList[nSNVType].apstRec[i],g_stSNVList[nSNVType].anUsed[i],sizeof(stRec), CompareSNVRecs); 
    }

    printf("GenerateSNVList[%d]...Done!!!\n",nSNVType);
}


int GetIndex(stRec *pstRec, int nTotalRec, int nOffset)
{
     int nStart=0, nLast=nTotalRec, nSeek; 

     do
     {
         nSeek = (nStart+nLast)/2;
         if (pstRec[nSeek].nOffset == nOffset){
             while(nSeek> 0 && pstRec[nSeek-1].nOffset == nOffset) nSeek--;
             return nSeek;   
         }
         else if (pstRec[nSeek].nOffset > nOffset) nLast = nSeek-1;
         else nStart = nSeek+1;
         
     }while(!(nStart > nLast));
     
     return -1;   
}


void GetSVarSeqVariant(int nChro, int nStart, int nEnd, const char* pcSeq, stSeqSNV &ostSeq)
{
     stSeqSNV astSeqSNV[3]; stRec *pstRec; ostSeq.nUsed=0; int nIdx;
  
     for (int nBp=nStart; nBp<nEnd+1; nBp++)
     {
         for (int nSNV=0; nSNV<3; nSNV++) {
              astSeqSNV[nSNV].nUsed=0; 
              nIdx = GetIndex(g_stSNVList[nSNV].apstRec[nChro],g_stSNVList[nSNV].anUsed[nChro],nBp);
              
              if (nIdx==-1) continue;

              do
              {                     
                   pstRec = g_stSNVList[nSNV].apstRec[nChro];
                      
                   if (nSNV == eSNP) astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[2]);
                   else if (nSNV == eINS) {
                           if (astSeqSNV[nSNV].nUsed > 0){astSeqSNV[nSNV].Add(',');}
                                astSeqSNV[nSNV].Add('+'); astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[0]);                       
                           }
                   else {astSeqSNV[nSNV].Add('-');astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[0]);}                 
                   nIdx++; 
              }while(!(nIdx > g_stSNVList[nSNV].anUsed[nChro]) && pstRec[nIdx].nOffset == nBp);
         }

         if (astSeqSNV[eSNP].nUsed==0 && astSeqSNV[eINS].nUsed==0 && astSeqSNV[eDEL].nUsed==0)
            {ostSeq.Add(pcSeq[nBp-nStart]); continue;}
         
         if (astSeqSNV[eSNP].nUsed==0 && astSeqSNV[eDEL].nUsed==0) ostSeq.Add(pcSeq[nBp-nStart]);  
 
         if (astSeqSNV[eSNP].nUsed > 0||astSeqSNV[eDEL].nUsed > 0) 
         {       
            ostSeq.Add('[');
            if (astSeqSNV[eSNP].nUsed > 0) {
                ostSeq.Add(pcSeq[nBp-nStart]); ostSeq.Add(astSeqSNV[eSNP].pcData,astSeqSNV[eSNP].nUsed); 
            }
    
            if (astSeqSNV[eDEL].nUsed > 0) {
                if (astSeqSNV[eSNP].nUsed > 0){ostSeq.Add(',');} 
                ostSeq.Add(astSeqSNV[eDEL].pcData,astSeqSNV[eDEL].nUsed);                   
            }
            ostSeq.Add(']');
         }

         if (astSeqSNV[eINS].nUsed > 0){
            ostSeq.Add('['); ostSeq.Add(astSeqSNV[eINS].pcData,astSeqSNV[eINS].nUsed); ostSeq.Add(']');
         }
    }    
}


void ProcessSeqRec()
{
    char acbuf[40960],*pChr=NULL,acTmp[4096],acOut1[40960],acOut2[40960]; 
    int nChro,nOffset,nLength,nUsed,nTotalSeq,nOutLen;       
    stSeqSNV ostSeq; stRec *pstRec;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; 
 
    while (!feof(g_pfSeq))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfSeq)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

        if (acbuf[0]==0||acbuf[0]=='#') continue;        

        if (acbuf[0]=='C'){
           fprintf(g_pfOut,"Chromosome1\tStart1\tEnd1\tLeft1\tRight1\tChromosome2\tStart2\tEnd2\tLeft2\tRight2\n",acbuf); 
           continue;
        }

        pChr = strtok(acbuf,"\t"); 

        if (pChr[0] == 'X') nChro = 22; else if (pChr[0] == 'Y') nChro = 23; else nChro = atoi(pChr)-1;          

        if (nChro < 22) {fprintf(g_pfOut,"%d\t",nChro+1); sprintf(acOut1,"%d\t%c",nChro+1,'\0');}
        else {fprintf(g_pfOut,"%c\t",pChr[0]); sprintf(acOut1,"%c\t%c",pChr[0],'\0');}

        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr); fprintf(g_pfOut,"%d\t",nOffset);        
        nOutLen = strlen(acOut1); sprintf(&acOut1[nOutLen],"%d\t%c",nOffset,'\0');

        pChr = strtok(NULL,"\t"); nLength = atoi(pChr); fprintf(g_pfOut,"%d\t",nOffset+nLength);
        nOutLen = strlen(acOut1); sprintf(&acOut1[nOutLen],"%d\t%c",nOffset+nLength,'\0');

        pChr = strtok(NULL,"\t");   
        GetSVarSeqVariant(nChro, nOffset-500, nOffset-1, pChr, ostSeq);   
        ostSeq.Add('\0'); fprintf(g_pfOut,"%s\t%s\t",ostSeq.pcData,pChr+500);             
        nOutLen = strlen(acOut1); sprintf(&acOut1[nOutLen],"%s\t%s\t%c",ostSeq.pcData,pChr+500,'\0');
         
        pChr = strtok(NULL,"\t"); nOffset+=nLength; nLength+=2; 
        GetSVarSeqVariant(nChro, nOffset, nOffset+500, pChr+nLength, ostSeq);
        ostSeq.Add('\0');
        
        memcpy(&acTmp[0],pChr,nLength); acTmp[nLength]=0;         
        fprintf(g_pfOut,"%s%s\t",acTmp,ostSeq.pcData);  
        nOutLen = strlen(acOut1); sprintf(&acOut1[nOutLen],"%s%s%c",acTmp,ostSeq.pcData,'\0');
         
        pChr = strtok(NULL,"\t");

        if (pChr[0] == 'X') nChro = 22; else if (pChr[0] == 'Y') nChro = 23; else nChro = atoi(pChr)-1;

        if (nChro < 22) {fprintf(g_pfOut,"%d\t",nChro+1); sprintf(acOut2,"%d\t%c",nChro+1,'\0');}
        else {fprintf(g_pfOut,"%c\t",pChr[0]); sprintf(acOut2,"%c\t%c",pChr[0],'\0');}

        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr); fprintf(g_pfOut,"%d\t",nOffset);
        nOutLen = strlen(acOut2); sprintf(&acOut2[nOutLen],"%d\t%c",nOffset,'\0');

        pChr = strtok(NULL,"\t"); nLength = atoi(pChr); fprintf(g_pfOut,"%d\t",nOffset+nLength);
        nOutLen = strlen(acOut2); sprintf(&acOut2[nOutLen],"%d\t%c",nOffset+nLength,'\0');  

        pChr = strtok(NULL,"\t");
        GetSVarSeqVariant(nChro, nOffset-500, nOffset-1, pChr, ostSeq);
        ostSeq.Add('\0'); fprintf(g_pfOut,"%s\t%s\t",ostSeq.pcData,pChr+500);
        nOutLen = strlen(acOut2); sprintf(&acOut2[nOutLen],"%s\t%s\t%c",ostSeq.pcData,pChr+500,'\0');

        pChr = strtok(NULL,"\t"); nOffset+=nLength; nLength+=2;
        GetSVarSeqVariant(nChro, nOffset, nOffset+500, pChr+nLength, ostSeq);
        ostSeq.Add('\0');

        memcpy(&acTmp[0],pChr,nLength); acTmp[nLength]=0;            
        fprintf(g_pfOut,"%s%s\n",acTmp,ostSeq.pcData);        
        nOutLen = strlen(acOut2); sprintf(&acOut2[nOutLen],"%s%s%c",acTmp,ostSeq.pcData,'\0'); 
 
        fprintf(g_pfOut,"%s\t%s\n",acOut2,acOut1);
//break;
    }
}


void CloseFiles()
{
    if (g_pfSeq) fclose(g_pfSeq); if (g_pfOut) fclose(g_pfOut);

}


int main(int argc, char *argv[])
{
    if (argc != 6){ printf("Invalid parameters\n");banner(argv); exit(9);}

    g_pfSeq = fopen(argv[1],"r"); if (!g_pfSeq){printf("Failed to open %s ...\n",argv[1]); goto ExitRtn;}

    for (int i=0; i<3; i++)
    {  
         g_stSNVList[i].pf = fopen(argv[i+2],"r"); 
         if (!g_stSNVList[i].pf){printf("Failed to open %s ...\n",argv[i+2]); goto ExitRtn;}
    }

    g_pfOut = fopen(argv[5],"w"); if (!g_pfOut){printf("Failed to open %s ...\n",argv[5]); goto ExitRtn;}
    
    for (int i=0; i<3; i++) GenerateSNVList(i); 

/*
    FILE *pfLog[3];

    pfLog[0]=fopen("Chk_snp.log","w"); pfLog[1]=fopen("Chk_ins.log","w"); pfLog[2]=fopen("Chk_del.log","w");

    stRec *pstRec; int nUsed, nChro;

    for (int i=0; i<3; i++){
        for (int nIdx=0; nIdx<24; nIdx++)
        {  
            nChro=nIdx+1;
            pstRec=g_stSNVList[i].apstRec[nIdx]; nUsed = g_stSNVList[i].anUsed[nIdx];  

            for (int j=0; j<nUsed; j++){
                 fprintf(pfLog[i],"%d\t%u\t%s\t%.4f\n",
                          nChro,pstRec[j].nOffset,pstRec[j].acAllele,pstRec[j].fCfdScore);  
            }
        }
    }    

    fclose(pfLog[0]); fclose(pfLog[1]); fclose(pfLog[2]);

*/   
    printf("Start Process Seq Recs ...\n"); ProcessSeqRec(); 
    printf("Done ...\n");     

ExitRtn:
    CloseFiles();     
} 
