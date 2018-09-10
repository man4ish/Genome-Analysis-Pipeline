/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on March 16, 2010, 11:04 AM
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
        fprintf(stdout,"Parameters not sufficient...\n");
        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename> <Output filename 1-greater than Min. SR> <Output filename 2-less than Min. SR > <Min. Supporting Reads> \n");
        exit(9);
    }
    FILE *pfIn, *pfOut, *pfFiltered;

    pfIn = fopen(argv[1],"r"); if (!pfIn) {fprintf(stdout,"Failed to open %s ...\n",argv[1]); exit(9);}
    pfOut = fopen(argv[2],"w");if (!pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[2]); exit(9);}
    pfFiltered = fopen(argv[3],"w");if (!pfFiltered) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); exit(9);}

    int nMinSR = atoi(argv[4]);

    char acbuf[40960],acOut[40960]; char *pChr=NULL; int nSR=0;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while (!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr=0;}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}

        if (acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0] ==0){
           continue;
        }

        strcpy(acOut,acbuf);
       
        pChr = strtok(acbuf,"\t"); //Chromosome
        pChr = strtok(NULL,"\t");  //Offset
        pChr = strtok(NULL,"\t");  //Neucleotide variant
        pChr = strtok(NULL,"\t");  //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  //Supporting Reads

        nSR = atoi(pChr);
        if (nSR < nMinSR) {
           fprintf(pfFiltered,"%s\n",acOut);
           continue;
         }
        fprintf(pfOut,"%s\n",acOut);
    }

    if (pfIn) fclose(pfIn);
    if (pfOut) fclose(pfOut);
    if (pfFiltered) fclose(pfFiltered);



    return (EXIT_SUCCESS);

}

