#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char** argv) {

    if (argc < 4)
    {
        fprintf(stdout,"Parameters not sufficient...\n");

        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename> <Tag :- SNP/INS/DEL> <Output Known dbsnp filename>\n");

        exit(9);
    }

    unsigned long ulTest1=0, ulTest2=0;
    FILE *pfIn, *pfOut; char acTag[4], acbuf[1024],*pChr=NULL;

    pfIn = fopen(argv[1],"r"); if (!pfIn) { fprintf(stdout,"Failed to opend %s ...\n", argv[1]); goto ExitRtn;}
    strcpy(acTag,argv[2]);
    pfOut = fopen(argv[3],"w"); if (!pfOut) {fprintf(stdout,"Failed to opend %s ...\n", argv[3]); goto ExitRtn;}  

    while (!feof(pfIn))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;
        
        pChr = strchr(acbuf,'\r'); if (pChr) *pChr='\0';
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr='\0';
        
        if (acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0]==0) continue;
        
        fprintf(pfOut,"%s\t%s\n",acbuf,acTag);
    }   

       
ExitRtn:
    if (pfIn) fclose(pfIn); if (pfOut) fclose(pfOut);    
}
