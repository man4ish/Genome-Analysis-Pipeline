#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *g_pfIn, *g_pfOut; int g_nGapSize=0;
char g_cType;

void ExtractRecs()
{
    char acbuf[4096],acOut[4096], acChro[3], acAllele[4], *pChr=NULL;
    int nChromosome, nOffset, nSR, nDiff; unsigned long ulTotalCnt;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"FilterRecs-Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}
       
        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t");//Chromosome
        pChr = strtok(NULL,"\t"); //SNP Info.
        pChr = strtok(NULL,"\t"); //Insert Info
        pChr = strtok(NULL,"\t"); //Insert Gap
        
        if (g_cType == 'I'){         
            if (atoi(pChr) != g_nGapSize) continue;
        }
        else if (g_cType == 'D'){        
            pChr = strtok(NULL,"\t"); //Delete Info
            //pChr = strtok(NULL,"\t"); //Delete Gap
            if (atoi(pChr) != g_nGapSize) continue;
        }      
        fprintf(g_pfOut,"%s\n",acOut);
    }
}


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n");
       fprintf(stdout,"\t<Input filename> <Gap Size> <Type:-I,D> <Output filename>\n\n");
       exit(9);
    }

    g_pfIn = fopen(argv[1],"r"); if (!g_pfIn) {fprintf(stdout,"Failed to open %s ...",argv[1]); goto ExitRtn;}
    g_nGapSize = atoi(argv[2]);  g_cType = argv[3][0];
    g_pfOut = fopen(argv[4],"w"); if (!g_pfOut) {fprintf(stdout,"Failed to open %s ...",argv[4]); goto ExitRtn;}

    ExtractRecs();    

ExitRtn:
  
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut); 
}
