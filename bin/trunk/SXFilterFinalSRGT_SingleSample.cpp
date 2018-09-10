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
        fprintf(stdout,"argv[0] = %s\n",argv[0]);
        fprintf(stdout,"argv[1] = %s\n",argv[1]);
        fprintf(stdout,"argv[2] = %s\n",argv[2]);
        fprintf(stdout,"argv[3] = %s\n",argv[3]);
 
        fprintf(stdout,"Parameters not sufficient...\n");
        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename> <Output filename> <Max. Supporting Reads>\n"); 
        exit(9);
    }
    FILE *pfIn, *pfOut;

    pfIn = fopen(argv[1],"r");
    pfOut = fopen(argv[2],"w");
    unsigned int unMaxSR = atoi(argv[3]);


    char acbuf[1024],acOut[1024]; char *pChr=NULL; short nTabCnt=0;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while (!feof(pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = 36.3=== %uM\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),pfIn)) break;

        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;

//        if (acbuf[0] == '\n'||acbuf[0] == '\r'){
//           fprintf(pfOut,"\n"); fprintf(pfOut,"\n"); continue;
//        }

        if (acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0] ==0 ){
           //fprintf(pfOut,"%s\n",acbuf);
           continue;
        }

        strcpy(acOut,acbuf); 
        pChr = strtok(acbuf,"\t"); //Chromosome
        pChr = strtok(NULL,"\t");  //Offset
        pChr = strtok(NULL,"\t");  //Nucleotide_Variant  
        pChr = strtok(NULL,"\t");  //Fwd_SNP_Reads
        pChr = strtok(NULL,"\t");  //Rvs_SNP_Reads
        pChr = strtok(NULL,"\t");  //Total_SNP_Reads

        if (atoi(pChr) > unMaxSR) continue;
        fprintf(pfOut,"%s\n",acOut);
    }

    if (pfIn) fclose(pfIn);
    if (pfOut) fclose(pfOut);



    return (EXIT_SUCCESS);

}

