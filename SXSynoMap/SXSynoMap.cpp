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
    printf("Synamatix SXSynoMap Copyright 2009 Synamatix Sdn Bhd\n");
    printf("Build %s %s\n",__DATE__,__TIME__);
    printf("Usage: %s <database name> <input file name> <output file name>\n\n", argv[0]);
}


int main(int argc, char** argv) {
    if (argc <4) {printf("Parameters required...\n");banner(argv); exit(9);}

    SXSynabase *pSXb = new SXSynabase(argv[1]);//("h_sapiens_genome_ncbi36v2-sb3.0.9");

    if (!pSXb) {fprintf(stderr,"Failed to create SXSynabase object...");exit(9);}

    FILE *pfSyn = fopen(argv[2],"r");//"SNP.syn","r"); 
    if (!pfSyn){fprintf(stderr,"Failed to open %s ...\n",argv[2]); exit(9);}

    FILE *pf=fopen(argv[3],"w");//"SNP.fna","w");

    if (!pf){fprintf(stderr,"Failed to open %s ...\n",argv[3]); exit(9);}
 
    char acbuf[256]; unsigned int nBufSize = sizeof(acbuf); 
    char *pChr; char *lpc[5]; int nIdx;  t_sindex sidx;  char acSeq[90];

    while(!feof(pfSyn))
    {
        if (!fgets(acbuf,nBufSize,pfSyn)) break;
        
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr=0;}
        nIdx=0; pChr = strtok(acbuf,"|");

        while (pChr)
        {
             lpc[nIdx++]= pChr;
             if (nIdx > 4) break;  pChr = strtok(NULL,"|");
        }
             
        //RecNo,GiNum,Offset,RefBase>VarBase 
        //fprintf(stderr,"%s\t%s\t%s\t%s\t%s\n",lpc[0],lpc[1],lpc[2],lpc[3],lpc[4]);                           
        if (pSXb->SXGetGiIdx(lpc[1],&sidx))
	{		
            if (pSXb->SXGetSequence(sidx,atol(lpc[2])-45,acSeq,89)> 0)	
            { 
               acSeq[44]=lpc[3][1];
               fprintf(pf,">%s|%s|%s|%s|%s\n%s\n",lpc[0],lpc[1],lpc[2],lpc[3],lpc[4],acSeq);
            }
        }
              
        memset(lpc,0,sizeof(lpc));        
    }

    if (pfSyn) fclose(pfSyn); if (pf) fclose(pf);
    if (pSXb) delete(pSXb);

    return (EXIT_SUCCESS);
}

