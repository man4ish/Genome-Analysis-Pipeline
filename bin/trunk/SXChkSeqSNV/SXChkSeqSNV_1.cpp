#include "SXChkSeqSNV.h"


static void banner(char *argv[])
{
    printf("Synamatix SXChkSeqSNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: \n", argv[0]);
    printf("\t<Seq filename> <SNP list filename> <INS list filename> <DEL list filename> <Output filename>\n\n");
}


void GenerateSNPList(int nSNVType)
{
    char acbuf[40960], *pChr=NULL,*pcLast=NULL; stRec *pstRec; 
    int nUsed, nIdx; 

    while (!feof(g_stSNVList[nSNVType].pf))
    {
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
}


void CheckSNVSeq(int nChro, int nStart, int nEnd, const char* pcSeq, stSeqSNV &ostSeq)
{
     int nTotalSeq; ostSeq.nUsed=0;
  
     for (int nBp=nStart; nBp < nEnd+1; nBp++)
     {
         for (int nSNV=0; nSNV<3; nSNV++) {
              pstRec = g_stSNVList[nSNV].apstRec[nChro];
              for (int i=0; i< g_stSNVList[nSNV].anUsed[nChro]; i++) {
                   if (pstRec[i].nOffset > nBp) break; if (pstRec[i].nOffset < nBp) continue;
                   if (pstRec[i].nOffset == nBp) {
                       /*if (astSeqSNV[nSNV].nUsed == astSeqSNV[nSNV].nTotal) {
                           astSeqSNV[nSNV].nTotal += SEQSNV_SIZE;
                           astSeqSNV[nSNV].pcData = (char*)realloc(astSeqSNV[nSNV].pcData,astSeqSNV[nSNV].nTotal);
                        }

                        memcpy(&(astSeqSNV[nSNV].pcData[astSeqSNV[nSNV].nUsed]),&(pstRec[i].acAllele[2]),1); astSeqSNV[nSNV].nUsed++;
                        */
                     pstRec[i].Add(&(pstRec[i].acAllele[2]),1);  
                   }
              }
              nTotalSeq += astSeqSNV[nSNV].nUsed;
         }
         /* 
         if (ostSeq.nUsed+nTotalSeq > ostSeq.nTotal){
             ostSeq.nTotal += (nTotalSeq+SEQSNV_SIZE);
             ostSeq.pcData = (char*)realloc(ostSeq.pcData,ostSeq.nTotal);
         }*/

         if (astSeqSNV[eSNP].nUsed==0 && astSeqSNV[eINS].nUsed==0 && astSeqSNV[eDEL].nUsed==0){
             ostSeq.pcData[ostSeq.nUsed]=pcSeq[nBp-nStart];
             ostSeq.nUsed++; continue;
         }

         ostSeq.pcData[ostSeq.nUsed]='['; ostSeq.nUsed++;

         if (astSeqSNV[eSNP].nUsed > 0){
             memcpy(&(ostSeq.pcData[ostSeq.nUsed]),&(astSeqSNV[eSNP].pcData),astSeqSNV[eSNP].nUsed);
             ostSeq.nUsed += astSeqSNV[eSNP].nUsed;
         }

         if (astSeqSNV[eINS].nUsed > 0){
             ostSeq.pcData[ostSeq.nUsed]='('; ostSeqSNV.nUsed++;
             memcpy(&(ostSeq.pcData[ostSeq.nUsed]),astSeqSNV[eINS].pcData,astSeqSNV[eINS].nUsed);
             ostSeq.nUsed += astSeqSNV[eINS].nUsed;
             ostSeq.pcData[ostSeq.nUsed]=')'; ostSeq.nUsed++;
         }

         if (astSeqSNV[eDEL].nUsed > 0){
             memcpy(&(ostSeq.pcData[ostSeq.nUsed]),astSeqSNV[eDEL].pcData,astSeqSNV[eDEL].nUsed);
             ostSeq.nUsed += astSeqSNV[eDEL].nUsed;
         }

         ostSeq.pcData[ostSeq.nUsed]=']'; ostSeq.nUsed++;
    }
    ostSeq.pcData[ostSeq.nUsed]='0'; 
}


void ProcessSeqRec()
{
    char acbuf[300],acOut[300],*pChr=NULL; int nChro,nOffset,nUsed,nTotalSeq; 
    stSeqSNV astSeqSNV[3],ostSeq; 
    stRec *pstRec;
 
    while (!feof(g_pfSeq))
    {
        if (!fgets(acbuf,sizeof(acbuf),g_pfSeq)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

        if (acbuf[0]==0||acbuf[0]=='#') continue;
        
        if (acbuf[0]=='C'){
            fprintf(g_pfOut,"%s\tLeft\tRight\n",acbuf);
            continue;
        }

        strcpy(acOut,acbuf); fprintf(g_pfOut,"%s\t",acOut);

        pChr = strtok(acbuf,"\t"); 

        if (pChr[0] == 'X') nChro = 22;
        else if (pChr[0] == 'Y') nChro = 23;
        else nChro = atoi(pChr);
 
        pChr = strtok(NULL,"\t"); nOffset = atoi(pChr);         
        pChr = strtok(NULL,"\t");   

        CheckSNVSeq(nChro, nOffset-100, nOffset+1, pChr, &ostSeq);   
        fprintf(g_pfOut,"%s\t",ostSeq.pcData);  
/*
        for (int nBp=nOffset-100; nBp < nOffset+1; nBp++)   
        {            
            for (int nSNV=0; nSNV<3; nSNV++) {                               
                pstRec = g_stSNVList[nSNV].apstRec[nChro];
                for (int i=0; i< g_stSNVList[nSNV].anUsed[nChro]; i++) {
                    if (pstRec[i].nOffset > nBp) break; if (pstRec[i].nOffset < nBp) continue; 
                    if (pstRec[i].nOffset == nBp) { 
                        if (astSeqSNV[nSNV].nUsed == astSeqSNV[nSNV].nTotal) {
                           astSeqSNV[nSNV].nTotal += SEQSNV_SIZE;
                           astSeqSNV[nSNV].pcData = (char*)realloc(astSeqSNV[nSNV].pcData,astSeqSNV[nSNV].nTotal);                          
                        }  

                        memcpy(&(astSeqSNV[nSNV].pcData[astSeqSNV[nSNV].nUsed]),&(pstRec[i].acAllele[2]),1); astSeqSNV[nSNV].nUsed++;
                    }
                }                  
                nTotalSeq += astSeqSNV[nSNV].nUsed; 
            }               

            if (ostSeq.nUsed+nTotalSeq > ostSeq.nTotal){
                ostSeq.nTotal += (nTotalSeq+SEQSNV_SIZE); 
                ostSeq.pcData = (char*)realloc(ostSeq.pcData,ostSeq.nTotal);                
            }
 
            if (astSeqSNV[eSNP].nUsed==0 && astSeqSNV[eINS].nUsed==0 && astSeqSNV[eDEL].nUsed==0){     
                ostSeq.pcData[ostSeq.nUsed]=pChr[ostSeq.nUsed];
                ostSeq.nUsed++; continue;   
            }

            ostSeq.pcData[ostSeq.nUsed]='['; ostSeq.nUsed++;

            if (astSeqSNV[eSNP].nUsed > 0){                
                memcpy(&(ostSeq.pcData[ostSeq.nUsed]),&(astSeqSNV[eSNP].pcData),astSeqSNV[eSNP].nUsed);            
                ostSeq.nUsed += astSeqSNV[eSNP].nUsed;                                             
            }

            if (astSeqSNV[eINS].nUsed > 0){
                ostSeq.pcData[ostSeq.nUsed]='('; ostSeqSNV.nUsed++; 
                memcpy(&(ostSeq.pcData[ostSeq.nUsed]),astSeqSNV[eINS].pcData,astSeqSNV[eINS].nUsed);
                ostSeqSNV.nUsed += astSeqSNV[eINS].nUsed;
                ostSeq.pcData[ostSeq.nUsed]=')'; ostSeq.nUsed++; 
            }

            if (astSeqSNV[eDEL].nUsed > 0){
                memcpy(&(ostSeq.pcData[ostSeq.nUsed]),astSeqSNV[eDEL].pcData,astSeqSNV[eDEL].nUsed);
                ostSeqSNV.nUsed += astSeqSNV[eDEL].nUsed;
            }

            ostSeq.pcData[ostSeq.nUsed]=']'; ostSeq.nUsed++;
        }
*/        
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
    
    for (int i=0; i<3; i++) GenerateSNPList(i); 


    FILE *pfLog=fopen("Chk.log","w"); 

    stRec *pstRec; int nUsed, nChro;

    for (int i=0; i<3; i++){
        for (int nIdx=0; nIdx<24; nIdx++)
        {  
            nChro=nIdx+1;
            pstRec=g_stSNVList[i].apstRec[nIdx]; nUsed = g_stSNVList[i].anUsed[nIdx];  

            for (int j=0; j<nUsed; j++){
                 fprintf(pfLog,"%d\t%u\t%s\t%.4f\n",
                          nChro,pstRec[j].nOffset,pstRec[j].acAllele,pstRec[j].fCfdScore);  
            }
        }
    }    

    fclose(pfLog);

    ProcessSeqRec(); 
    

ExitRtn:
    CloseFiles();     
} 
