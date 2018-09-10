#include "SXMergeMFiles.h"

M_List g_MList;


FILE *g_pfIn1=NULL, *g_pfIn2=NULL, *g_pfOut=NULL;

static void banner(char *argv[])
{
    printf("Synamatix SXMergeMFiles Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename1- ANU final list> <Input filename2 - ANU aqs list> <Output filename>\n\n");
}


void GenerateInFileList()
{
    fprintf(stdout,"Generating InFile List...\n");

    char acbuf[4096],*pChr=NULL; stKey *pKey=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn1))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn1)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C') continue;

        pChr = strtok(acbuf,"\t"); //Chromosome

        pChr = strtok(NULL,"\t");  //Offset
        pKey=new stKey; pKey->nOffset = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Allele
        pKey->pcAllele = new char[strlen(pChr)+1];
        strcpy(pKey->pcAllele,pChr);

        g_MList[pKey]=1;            
   }   
}


void ProcessRecs()
{
     fprintf(stdout,"ProcessRecs...\n");
     M_List::iterator itr;

     char acbuf[4096], acOut[4096]; char *pChr=NULL; stKey *pKey=new stKey;
     unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

     while (!feof(g_pfIn2))
     {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %u M\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn2)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'||acbuf[0] == 'C'){
           fprintf(g_pfOut,"%s\n",acbuf); continue;}

        strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t"); //Chromosome  
        pChr = strtok(NULL,"\t");  //Offset
        pKey->nOffset = atoi(pChr);

        pChr = strtok(NULL,"\t");  //Allele
        pKey->pcAllele = new char[strlen(pChr)+1];
        strcpy(pKey->pcAllele,pChr);

        itr = g_MList.find(pKey);

        if (itr != g_MList.end()){
            fprintf(g_pfOut,"%s\n",acOut);
        }  
     }
}


void CloseFiles()
{
    if (g_pfIn1) fclose(g_pfIn1); if (g_pfIn2) fclose(g_pfIn2); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv){

    if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfIn1 = fopen(argv[1],"r");
    if (!g_pfIn1){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfIn2 = fopen(argv[2],"r");
    if (!g_pfIn2){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    GenerateInFileList(); ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);

}


