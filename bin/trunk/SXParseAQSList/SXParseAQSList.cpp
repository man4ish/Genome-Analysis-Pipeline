/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on March 8, 2010, 11:56 AM
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 4)
    {
        fprintf(stdout,"argv[0]=%s\n",argv[0]);
        fprintf(stdout,"argv[1]=%s\n",argv[1]);
        fprintf(stdout,"argv[2]=%s\n",argv[2]);
        fprintf(stdout,"argv[3]=%s\n",argv[3]);
    
        fprintf(stdout,"Parameters not sufficient...\n");

        fprintf(stdout,"Usage: %s\n",argv[0]);  
        fprintf(stdout,"\t<Input filename> <Output Known dbsnp filename> <Output Novel filename>\n");
    
        exit(9);
    }
    FILE *pfIn, *pfOutdbSnp, *pfOutdbNoval;

    pfIn = fopen(argv[1],"r"); if (!pfIn) {fprintf(stdout,"Failed to opend %s ...\n", argv[1]); exit(9);}
    pfOutdbSnp = fopen(argv[2],"w");
    pfOutdbNoval = fopen(argv[3],"w");

    char acbuf[4096]; char *pChr=NULL; short nTabCnt=0;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while (!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %lu x 0.1M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;

        if (acbuf[0] == '\n'){
           fprintf(pfOutdbSnp,"\n"); fprintf(pfOutdbNoval,"\n"); continue;
        }

        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '#'||acbuf[0] == 'C'){
           fprintf(pfOutdbSnp,"%s\n",acbuf);
           fprintf(pfOutdbNoval,"%s\n",acbuf);continue;
        }

        pChr = strchr(acbuf,'\t');
        nTabCnt=0;
        do
        {
            pChr = strchr(pChr+1,'\t');
            nTabCnt++;

        } while (nTabCnt < 9);

        if (*(pChr+1)== '-')
           fprintf(pfOutdbNoval,"%s\n",acbuf);
        else
            fprintf(pfOutdbSnp,"%s\n",acbuf);
    }

    if (pfIn) fclose(pfIn);
    if (pfOutdbSnp) fclose(pfOutdbSnp);
    if (pfOutdbNoval) fclose(pfOutdbNoval);



    return (EXIT_SUCCESS);
}

