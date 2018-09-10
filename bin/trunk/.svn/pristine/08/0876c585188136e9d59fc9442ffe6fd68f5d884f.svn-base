#include "SXFilterSNV.h"


static void banner(char *argv[])
{
    printf("Synamatix SXFilterSNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input SNV filename> <Input Clustered filename with data to be filtered>"
           " <Output SNV filename> <Output SNV filename which will stored data filtered>\n\n");
}


void CloseFiles()
{
    if (g_fd > 0) close(g_fd);                   
    if (g_pfIn) fclose(g_pfIn); if (g_pfFilterIn) fclose(g_pfFilterIn);
    if (g_pfOut) fclose(g_pfOut); if (g_pfFilterOut) fclose(g_pfFilterOut);
}


void GenerateFilteredList()
{
    char acbuf[1024],*pChr=NULL; stKey *pstKey=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfFilterIn))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %lu x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfFilterIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pstKey = new stKey; 
        
        pChr = strtok(acbuf,"\t"); pstKey->nChromo = atoi(pChr);
        pChr = strtok(NULL,"\t");  pstKey->nOffset = atoi(pChr);        
        
        g_Itr = g_FilterList.find(pstKey);

        if (g_Itr == g_FilterList.end()) g_FilterList[pstKey]=1;        
        else delete[] pstKey;                 
   } 
}


inline unsigned int ComputeOffset(unsigned char acVal[4])
{
    unsigned int unOffset= acVal[3];unOffset<<=8;
    unOffset|= acVal[2];unOffset<<=8;
    unOffset|= acVal[1];unOffset<<=8;
    unOffset|= acVal[0];

    return unOffset;
}


void ProcessRecs()
{
    struct stat st; fstat(g_fd, &st);
    size_t file_size = (size_t)st.st_size; 
    unsigned long ulTotalRecs = file_size/sizeof(stSNV);

    stSNV* pSNV=(stSNV*)mmap(NULL,file_size,PROT_READ,MAP_PRIVATE,g_fd,0);
    
    if (pSNV == (void *)-1) {
        fprintf(stderr, "Error: mmap failed: %s\n", strerror(errno));return;
    }

    stKey *pstKey=new stKey;
 
    for (unsigned long i=0; i<ulTotalRecs; i++)
    {
        pstKey->nChromo = pSNV[i].nChromosome;
        pstKey->nOffset = ComputeOffset(pSNV[i].acOffset);

        g_Itr = g_FilterList.find(pstKey);

        if (g_Itr == g_FilterList.end()) fwrite(&pSNV[i],sizeof(stSNV),1,g_pfOut);
        else fwrite(&pSNV[i],sizeof(stSNV),1,g_pfFilterOut); 
    }

    int rc = munmap((char*)pSNV, file_size);
    if (rc == -1) fprintf(stderr, "munmap failed...\n");     

/*
    char acbuf[1024],acOut[1024],*pChr=NULL; stKey *pstKey=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn))
    {
        if (++ulTest1 == 100000)
        {
            fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pstKey = new stKey; strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t"); pstKey->cChromo = pChr[0];
        pChr = strtok(NULL,"\t");  pstKey->nOffset = atoi(pChr);

        g_Itr = g_FilterList.find(pstKey);

        if (g_Itr == g_FilterList.end()) fprintf(g_pfOut,"%s\n",acOut);        
        else fprintf(g_pfFilterOut,"%s\n",acOut);

        delete[] pstKey;        
   }    
*/

 delete[] pstKey;
}


int main(int argc, char** argv) {

    if (argc < 5){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_fd = open(argv[1],O_RDONLY, S_IRUSR);

    if (g_fd == -1){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfFilterIn = fopen(argv[2],"r");
    if (!g_pfFilterIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[3],"wb");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    g_pfFilterOut = fopen(argv[4],"wb");
    if (!g_pfFilterOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[4]); goto Err_Exit;
    }

    GenerateFilteredList(); ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}
