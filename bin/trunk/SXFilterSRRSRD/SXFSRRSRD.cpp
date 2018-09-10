#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
using namespace std;

int main (int argc, char* argv[])
{
  if (argc < 6){
     fprintf(stdout,"argv[0]=%s\n",argv[0]);
     fprintf(stdout,"argv[1]=%s\n",argv[1]);
     fprintf(stdout,"argv[2]=%s\n",argv[2]);
     fprintf(stdout,"argv[3]=%s\n",argv[3]);
     fprintf(stdout,"argv[4]=%s\n",argv[4]);
 
     fprintf(stdout,"Not sufficient parameters...\n");  
     fprintf(stdout,"Usage: %s\n",argv[0]);
     fprintf(stdout,"\t<Input filename><Min. Supporting Reads> <Min. Read Strength> "
                    "<Min. Read Density> <Output filename>\n");
     exit(9);    
  }

  FILE *pf = fopen(argv[1],"r"); if (!pf) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); exit(9);}
  int nSR = atoi(argv[2]);
  int nRS = atoi(argv[3]);
  int nRD = atoi(argv[4]);
  FILE *pfOut = fopen(argv[5],"w"); if (!pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[5]); exit(9);}
    
  char acbuf[4096],acOut[4096], *pChr=NULL; 
  int nLen=sizeof(acbuf), nTSR, nTRD, nTRS;
  unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
  char acMsg[1024];

  while (!feof(pf))
  {
     if (++ulTest1 == 1000000)
     {
         ulTest2++;
         
         fprintf(stderr,"Total Recs Processed = %u M\n", ulTest2);
         //sprintf(acMsg,"Total Recs Processed = %u M\n", ulTest2);
         //cout << acMsg;
         ulTest1=0;
     }
    
     if (!fgets(acbuf,nLen,pf)) break;        

     pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
     pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

     if (acbuf[0]=='C'||acbuf[0]=='#'||acbuf[0]==0) continue;

     strcpy(acOut,acbuf);
     
     pChr = strtok(acbuf,"\t");pChr = strtok(NULL,"\t");pChr = strtok(NULL,"\t");
     pChr = strtok(NULL,"\t");pChr = strtok(NULL,"\t");pChr = strtok(NULL,"\t");
     nTSR = atoi(pChr);

     pChr = strtok(NULL,"\t");pChr = strtok(NULL,"\t");pChr = strtok(NULL,"\t");
     nTRD = atoi(pChr);

     if (nTRD < 1) continue;

     nTRS = (int)((nTSR*100/nTRD)+0.5);  

     if (nTSR >= nSR && nTRS >= nRS && nTRD >= nRD){
         fprintf(pfOut,"%s\n",acOut);     
     }      
  }
 
  if (pf) fclose(pf); if (pfOut) fclose(pfOut);
}
