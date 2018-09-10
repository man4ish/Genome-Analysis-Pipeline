/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on April 26, 2010, 3:29 PM
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
 * 
 */

int main(int argc, char** argv) {

    if (argc < 5){
        fprintf(stdout,"argv[0]=%s\n",argv[0]);
        fprintf(stdout,"argv[1]=%s\n",argv[1]);
        fprintf(stdout,"argv[2]=%s\n",argv[2]);
        fprintf(stdout,"argv[3]=%s\n",argv[3]);
        fprintf(stdout,"argv[4]=%s\n",argv[4]);

        fprintf(stdout,"Not sufficient parameters...\n");
        fprintf(stdout,"Usage: %s\n", argv[0]);
        fprintf(stdout,"\t<Input filename><Min. Supporting Reads><Min. Read Strength><Min. Read Density>\n");
        exit(9);
    }

    FILE *pfIn = fopen(argv[1],"r");
    if (!pfIn) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); exit(9);}

    int nSR = atoi(argv[2]);
    int nRS = atoi(argv[3]);
    int nRD = atoi(argv[4]);

    int nInSR=0, nInRS=0, nInRD=0;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    char acbuf[5000], acOut[5000]; char *pChr=NULL;

    //FILE *pfTest=fopen("test.lst","w");

    while (!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;
        
        if (acbuf[0]=='#'||acbuf[0]==0||acbuf[0]=='C') continue; 

        strcpy(acOut,acbuf);

        pChr=strtok(acbuf,"\t"); //Chromosome
        pChr=strtok(NULL,"\t"); //Offset
        pChr=strtok(NULL,"\t"); //Neucleotide variant
        pChr=strtok(NULL,"\t"); //Supporting Reads
        nInSR = atoi(pChr);
        pChr=strtok(NULL,"\t"); //Read Density
        nInRD = atoi(pChr); 

        if (nInRD==0) continue;
//int strength = (int)((float)supp/(float)density * 100.0 + 0.5);
        nInRS = (int)((float)nInSR/(float)nInRD*100.0+0.5);

        if (nInRS > 100 || nInRS < 0) nInRS = 101;
                   

        //fprintf(pfTest,"nInSR = %d\tnInRD = %d\tnInRS=%d\n",nInSR,nInRD,nInRS);  
                           
        if (nInSR >= nSR && nInRS >= nRS && nInRD >= nRD)
            fprintf(stdout,"%s\n",acOut);
    }

    if (pfIn) fclose(pfIn);
    //if (pfTest) fclose(pfTest);

    return (EXIT_SUCCESS);
}

