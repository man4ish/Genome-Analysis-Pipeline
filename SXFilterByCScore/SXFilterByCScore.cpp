#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

FILE *g_pfIn=NULL, *g_pfOut=NULL; float g_fFilterVal;


static void banner(char *argv[])
{
    printf("Synamatix SXAppendScore Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Filter Value - decimal> <Output filename>\n\n");
}


void ProcessRecs()
{
    char acbuf[40960],acOut[40960],*pChr=NULL, *pcLast=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; float fCScore;

    while(!feof(g_pfIn))
    {        
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %lu x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t");  

        while(pChr){
            pcLast = pChr; pChr = strtok(NULL,"\t");            
        } 

        fCScore = atof(pcLast);

        if ( fCScore < g_fFilterVal)  continue;

        fprintf(g_pfOut,"%s\n",acOut);
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv) { 

    if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_fFilterVal = atof(argv[2]);

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

