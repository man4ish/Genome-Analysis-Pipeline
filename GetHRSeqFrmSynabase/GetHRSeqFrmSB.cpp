/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on December 8, 2009, 10:48 AM
 */

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "/opt/synamatix/synabase/synabase/include/SXSynabase.h"
/*
 * 
 */

static void banner(char *argv[])
{
    printf("Synamatix GetHRSeqFrmSB Copyright 2009 Synamatix Sdn Bhd\n");
    printf("Build %s %s\n",__DATE__,__TIME__);
    printf("Usage: %s <database name>  <chromosome> <offset> <size>\n\n", argv[0]);
}


unsigned int g_unGINums_G36[24] = {89161185,89161199,89161205,89161207,51511721,89161210,
                                   89161213,51511724,89161216,89161187,51511727,89161190,
                                   51511729,51511730,51511731,51511732,51511734,51511735,
                                   42406306,51511747,51511750,89161203,89161218,89161220};


unsigned int g_unGINums_G37[24] = {224589800,224589811,224589815,224589816,224589817,224589818,
                                   224589819,224589820,224589821,224589801,224589802,224589803,
                                   224589804,224589805,224589806,224589807,224589808,224589809,
                                   224589810,224589812,224589813,224589814,224589822,224589823};


int main(int argc, char** argv) {

    if (argc <4) {printf("Parameters required...\n");banner(argv); exit(9);}
           
    unsigned int *punGINums=NULL; char *pChr=NULL;

    SXSynabase *pSXb = new SXSynabase(argv[1]); //"h_sapiens_genome_ncbi36v2-sb3.0.9");

    if (!pSXb) {fprintf(stdout,"Failed to create SXSynabase object...");exit(9);}

    pChr = strchr(argv[1],'7'); punGINums=(!pChr)?g_unGINums_G36:g_unGINums_G37;
            
    int nChromosome = atoi(argv[2]);    
    int nOffset = atoi(argv[3])-1;
    int nSize = atoi(argv[4]);
    int nGINum = punGINums[nChromosome-1]; 
     
    char acGINum[10]; sprintf(acGINum,"%d",nGINum);  t_sindex sidx;  char acSeq[nSize+1]; 

    if (pSXb->SXGetGiIdx(acGINum,&sidx))
    {		
       if (pSXb->SXGetSequence(sidx,nOffset,acSeq,nSize)> 0)	
       { 
          fprintf(stdout,"%s\n",acSeq);
       }
    }
    else
   {
        fprintf(stdout,"SXGetGiIdx return false. GiNumm = %s\n",acGINum);
   }
           
    return (EXIT_SUCCESS);
}

