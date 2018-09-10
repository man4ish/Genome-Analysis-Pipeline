#include <ctype.h>
#include "SXChkSeqSNV_SmallCap.h"


static void banner(char *argv[])
{
    printf("Synamatix SXChkSeqSNV_SmallCap Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: \n", argv[0]);
    printf("\t<Seq filename> <SNP list filename> <INS list filename> <DEL list filename> <Output filename> <Size> <Type:-S/I/D>\n\n");
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
     int nStart, nLast, nSeek;
     nStart=0; nLast = nTotalRec; 

/*
if (nOffset == 232457857){
   printf("Last=%d\n",nLast);
}
*/     

     do
     {
         nSeek = (nStart+nLast)/2;
         if (pstRec[nSeek].nOffset == nOffset){
             while(nSeek> 0 && pstRec[nSeek-1].nOffset == nOffset) nSeek--;
             return nSeek;   
         }
         else if (pstRec[nSeek].nOffset > nOffset){
/*
if (nOffset == 232457857){
   printf("Seek=%d\tOffset=%d\n",nSeek,pstRec[nSeek].nOffset);
}
*/            nLast = nSeek-1;}
         else{
/*
if (nOffset == 232457857){
   printf("Seek=%d\tOffset=%d\n",nSeek,pstRec[nSeek].nOffset);
}             
*/
            nStart = nSeek+1;}
         
     }while(!(nStart > nLast));
/*
if (nOffset == 23023232){
   printf("Start=%d\tLast=%d\tSeek=%d\n",nStart,nLast,nSeek);
}
*/
     return -1;   
}


void GetSNVSeqVariant(int nChro, int nStart, int nEnd, const char* pcSeq, stSeqSNV &ostSeq, const char cDir)
{
     stSeqSNV astSeqSNV[3]; stRec *pstRec; ostSeq.nUsed=0; int nIdx;
  
     for (int nBp=nStart; nBp<nEnd+1; nBp++)
     {
         for (int nSNV=0; nSNV<3; nSNV++) {
              astSeqSNV[nSNV].nUsed=0; 
              nIdx = GetIndex(g_stSNVList[nSNV].apstRec[nChro],g_stSNVList[nSNV].anUsed[nChro],nBp);

              if (nIdx==-1) continue;

              do {                     
                   pstRec = g_stSNVList[nSNV].apstRec[nChro];
                      
                   if (nSNV == eSNP) astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[2]);                                        
                   else if (nSNV == eINS) {
                        if (astSeqSNV[nSNV].nUsed > 0)  astSeqSNV[nSNV].Add(','); 
                        astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[0]);                                  
                   }
                   else astSeqSNV[nSNV].Add(pstRec[nIdx].acAllele[0]);                 
                   nIdx++; 
              } while(!(nIdx > g_stSNVList[nSNV].anUsed[nChro]) && pstRec[nIdx].nOffset == nBp);
         } //end for

         if (astSeqSNV[eSNP].nUsed == 0 && astSeqSNV[eDEL].nUsed == 0 && astSeqSNV[eINS].nUsed == 0)
         {    ostSeq.Add(pcSeq[nBp-nStart]); continue;  }

         if (((cDir=='L') && (nBp==nEnd))||((cDir=='R') && (nBp==nStart))) 
         { 
              if (g_cType == 'I'){
                  if (cDir=='L'){
                      if (astSeqSNV[eSNP].nUsed > 0||astSeqSNV[eDEL].nUsed > 0) ostSeq.Add(tolower(pcSeq[nBp-nStart]));    
                      else ostSeq.Add(pcSeq[nBp-nStart]);
                  }

                  ostSeq.Add('+');

                  if (astSeqSNV[eINS].nUsed > 1) 
                      ostSeq.Add(tolower(astSeqSNV[eINS].pcData[0]));               
                  else ostSeq.Add(astSeqSNV[eINS].pcData[0]);

                  ostSeq.Add('+');   
              }
              else {
                  if (g_cType == 'S'){  
                      if (astSeqSNV[eINS].nUsed > 0 || astSeqSNV[eDEL].nUsed > 0 || astSeqSNV[eSNP].nUsed > 1)                        
                          ostSeq.Add(tolower(pcSeq[nBp-nStart]));
                      else ostSeq.Add(pcSeq[nBp-nStart]);                    
                  }
                  else if (g_cType == 'D'){ 
                      if (astSeqSNV[eSNP].nUsed > 0 || astSeqSNV[eINS].nUsed > 0 || astSeqSNV[eDEL].nUsed > 1) 
                          ostSeq.Add(tolower(pcSeq[nBp-nStart]));
                      else ostSeq.Add(pcSeq[nBp-nStart]);                
                  }

                  if (cDir=='R' && astSeqSNV[eINS].nUsed > 0){
                      ostSeq.Add('+');
                      for (int i=0; i<astSeqSNV[eINS].nUsed; i++) { ostSeq.Add(astSeqSNV[eINS].pcData[i]); }
                      ostSeq.Add('+');
                  }
              }
              continue;
         }
         
         if (astSeqSNV[eSNP].nUsed > 0||astSeqSNV[eDEL].nUsed > 0)  ostSeq.Add(tolower(pcSeq[nBp-nStart])); 

         if (astSeqSNV[eINS].nUsed > 0) {
             ostSeq.Add('+'); for (int i=0; i<astSeqSNV[eINS].nUsed; i++) ostSeq.Add(astSeqSNV[eINS].pcData[i]); ostSeq.Add('+'); 
         }
    }    
}


void ProcessSeqRec()
{
    char acbuf[40960],*pChr=NULL; int nChro,nOffset,nUsed,nTotalSeq; 
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
            fprintf(g_pfOut,"Chromosome\tStart\tEnd\tLeft\tRight\n");
            continue;
        }

        pChr = strtok(acbuf,"\t"); 

        if (pChr[0] == 'X') nChro = 22; else if (pChr[0] == 'Y') nChro = 23; else nChro = atoi(pChr)-1;          
        if (nChro < 22) fprintf(g_pfOut,"%d\t",nChro+1); else fprintf(g_pfOut,"%c\t",pChr[0]);

//if  (nChro!=1) continue;
 
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr); fprintf(g_pfOut,"%d\t%d\t",nOffset,nOffset);        

//if (nOffset == 232457857){
//           
//}

        pChr = strtok(NULL,"\t");   
        GetSNVSeqVariant(nChro, nOffset-g_nSize, nOffset, pChr, ostSeq,'L');   
        ostSeq.Add('\0'); fprintf(g_pfOut,"%s\t",ostSeq.pcData);  
         
        pChr = strtok(NULL,"\t");
        GetSNVSeqVariant(nChro, nOffset, nOffset+g_nSize, pChr, ostSeq,'R');
        ostSeq.Add('\0');
 
        fprintf(g_pfOut,"%s\n",ostSeq.pcData);  
    }
}


void CloseFiles()
{
    if (g_pfSeq) fclose(g_pfSeq); if (g_pfOut) fclose(g_pfOut);

}


int main(int argc, char *argv[])
{
    if (argc != 8){ printf("Invalid parameters\n");banner(argv); exit(9);}

    g_pfSeq = fopen(argv[1],"r"); if (!g_pfSeq){printf("Failed to open %s ...\n",argv[1]); goto ExitRtn;}

    for (int i=0; i<3; i++)
    {  
         g_stSNVList[i].pf = fopen(argv[i+2],"r"); 
         if (!g_stSNVList[i].pf){printf("Failed to open %s ...\n",argv[i+2]); goto ExitRtn;}
    }

    g_pfOut = fopen(argv[5],"w"); if (!g_pfOut){printf("Failed to open %s ...\n",argv[5]); goto ExitRtn;}
    g_nSize = atoi(argv[6]); g_cType = argv[7][0];
    
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
    printf("Start Processing Seq Recs ...\n"); ProcessSeqRec(); 
    printf("Done ...\n");     

ExitRtn:
    CloseFiles();     
} 
