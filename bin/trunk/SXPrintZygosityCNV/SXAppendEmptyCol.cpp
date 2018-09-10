/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on June 1, 2010, 11:06 AM
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
/*
 * 
 */


static void banner(char *argv[])
{
    printf("Synamatix SXAppendEmptyCol Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Output filename>\n\n");
}


int main(int argc, char** argv) {

    FILE *pfIn=NULL, *pfOut=NULL;
    char acbuf[1024], *pChr=NULL; int nLen = sizeof(acbuf);

    if (argc < 3){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    pfIn = fopen(argv[1],"r");
    if (!pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    pfOut = fopen(argv[2],"w");
    if (!pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    while(!feof(pfIn))
    {
        if (!fgets(acbuf,nLen,pfIn)){break;}        

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

        if (acbuf[0]=='#' || acbuf[0] == 'C'|| acbuf[0]==0) continue;

        fprintf(pfOut,"%s\t-\n",acbuf);
    }
    
    return (EXIT_SUCCESS);    

Err_Exit:
    if (pfIn) fclose(pfIn); if (pfOut) fclose(pfOut);   
    banner(argv); exit(9);
}

