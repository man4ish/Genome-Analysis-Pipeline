#include "SXAppendCytoband.h"

static void banner(char *argv[])
{
     printf("Synamatix SXAppendCytoband Copyright 2011 Synamatix Sdn Bhd\n");
     printf("Built %s %s\n", __DATE__, __TIME__);
     printf("Usage: \n", argv[0]);
     printf("\t<Input filename> <Input Cytoband Ref. filename> <Output filename>\n\n");
}


bool GenerateCytobandsMap()
{
    char acbuf[40960],acInfo[11],*pChr=NULL; 
    int nChro, nStart, nStop;    
      
    while (!feof(g_pfCytobands))
    {        
        if (!fgets(acbuf,sizeof(acbuf),g_pfCytobands)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0; 

        pChr = strtok(acbuf,"\t");
        
        nChro = atoi(pChr+3); 
        if (nChro == 0) {
           if (pChr[3]=='X') nChro=23;
           else if (pChr[3]=='Y') nChro=24;
        }

        sprintf(acInfo,"%d",nChro);

        pChr = strtok(NULL,"\t");
        nStart = atoi(pChr);

        pChr = strtok(NULL,"\t");
        nStop = atoi(pChr);

        pChr = strtok(NULL,"\t");
        strcat(acInfo,pChr);
   
        g_Cytobands_Itr = g_Cytobands_Map.find(nChro);

        if (g_Cytobands_Itr == g_Cytobands_Map.end())
            g_Cytobands_Itr = g_Cytobands_Map.insert(std::pair<int,Cytoband_Chromosome*>(nChro,new Cytoband_Chromosome())).first;  

        g_Cytobands_Itr->second->insertCytoband(nStart,nStop,acInfo);
    }

    for (g_Cytobands_Itr = g_Cytobands_Map.begin(); g_Cytobands_Itr != g_Cytobands_Map.end(); g_Cytobands_Itr++) 
    {
         g_Cytobands_Itr->second->sortCytobands();
    }
}


void OutputCytobandInfo(int nChro, int nOffset, FILE *pf)
{
    char *pcInfo[20];
    g_Cytobands_Itr = g_Cytobands_Map.find(nChro);

    unsigned unCnt = g_Cytobands_Itr->second->search(nOffset, pcInfo,20);
    
    if (unCnt > 0) 
    {
        fprintf(pf,"%s", pcInfo[0]); 
        for (unsigned i = 1; i < unCnt; i++) fprintf(pf,",%s", pcInfo[i]);                    
    }
    else 
    {
        fprintf(pf,"-");
    }
}


void ProcessRecs()
{
    char acbuf[40960],acOut[40960],*pChr=NULL; 
    int nChro, nOffset;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;      

    while (!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stdout,"Total Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf, sizeof(acbuf),g_pfIn)) break;

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0]==0||acbuf[0]=='#') { fprintf(g_pfOut,"%s\n",acbuf); continue; }

        if (acbuf[0]=='C'){ fprintf(g_pfOut,"%s\tCytoband\n",acbuf); continue; }

        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t"); //Chromosome

        nChro = atoi(pChr);

        if (nChro == 0){ 
            if (pChr[0]=='X') nChro = 23; 
            else if (pChr[0]=='Y') nChro = 24;
        }

        pChr = strtok(NULL,"\t");  //GiNumber
        pChr = strtok(NULL,"\t");  //Offset
        nOffset = atoi(pChr);

        fprintf(g_pfOut,"%s\t",acOut); 
        OutputCytobandInfo(nChro,nOffset,g_pfOut);
        fprintf(g_pfOut,"\n");
    }
}


int main(int argc, char *argv[])
{
    if (argc != 4){ printf("Invalid parameters\n");banner(argv); exit(9);}

    g_pfIn = fopen(argv[1],"r");
    if (!g_pfIn) { fprintf(stdout,"Failed to open %s ....",argv[1]); goto ExitRtn;}

    g_pfCytobands = fopen(argv[2],"r");
    if (!g_pfCytobands) { fprintf(stdout,"Failed to open %s ....",argv[2]); goto ExitRtn;}
    
    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut) { fprintf(stdout,"Failed to open %s ....",argv[3]); goto ExitRtn;}

    GenerateCytobandsMap(); ProcessRecs();  

ExitRtn:
    if (g_pfIn) fclose(g_pfIn);  if (g_pfCytobands) fclose(g_pfCytobands);
    if (g_pfOut) fclose(g_pfOut); g_pfOut=NULL;
}
