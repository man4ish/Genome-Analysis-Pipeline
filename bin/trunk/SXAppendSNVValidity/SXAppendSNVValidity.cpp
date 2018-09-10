#include <stdio.h>
#include <stdlib.h>
#include <map>

typedef struct _stKey
{
    char *pcRsid;

    ~_stKey()
    {
         if (pcRsid) delete[] pcRsid;
    }
} stKey;


typedef struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
         return strcmp(p1->pcRsid,p2->pcRsid) < 0;
    }	
};

typedef struct _stData
{
    char *pcValidity;

    ~_stData()
    {
        if (pcValidity) delete[] pcValidity;   
    }
} stData;


typedef std::map<stKey*,stData*,_compareKey> VALIDITY_LIST;

VALIDITY_LIST g_Validity_List[24]; VALIDITY_LIST::iterator itr=NULL;
FILE *g_pfValidity=NULL, *g_pfIn=NULL, *g_pfOut=NULL;
unsigned long ulCnt=0;

FILE *pf = fopen("Check.lst","w");

void GenerateValidityList()
{
    fprintf(stdout,"Generating ValidityInList...\n");

    char acbuf[81920],acOut[81920],*pChr=NULL; stKey *pstKey; stData *pstData; VALIDITY_LIST::iterator itr;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; int nChroIdx=0;

    while(!feof(g_pfValidity))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stdout,"Total Recs Processed = %uM\n", ++ulTest2);
            ulTest1=0;
        }
       
        if (!fgets(acbuf,sizeof(acbuf),g_pfValidity)){break;}
        strcpy(acOut,acbuf);

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pChr = strtok(acbuf,"\t"); //Chromosome        

        if (strlen(pChr) > 5) continue;

        //if (ulTest2 == 24) fprintf(pf,"%s\n",acOut);
        //if (ulTest2 == 25) return;

        if (isdigit(pChr[3])) nChroIdx = atoi(pChr+3)-1;            
        else if (pChr[3] == 'X') nChroIdx = 22;
        else if (pChr[3] == 'Y') nChroIdx = 23;

        pChr = strtok(NULL,"\t");  //Start
        pChr = strtok(NULL,"\t");  //Stop
        pChr = strtok(NULL,"\t");  //Validity

        pstData = new stData;
        pstData->pcValidity = new char[strlen(pChr)+1];
        strcpy(pstData->pcValidity,pChr);

        pChr = strtok(NULL,"\t");  //Key
        pstKey = new stKey;
        pstKey->pcRsid = new char[strlen(pChr)+1];
        strcpy(pstKey->pcRsid,pChr); 

        itr = g_Validity_List[nChroIdx].find(pstKey);

        if (itr == g_Validity_List[nChroIdx].end()) {g_Validity_List[nChroIdx][pstKey]=pstData;ulCnt++;}

        //else fprintf(stdout,"%s\n",pstKey->pcRsid);
    }
    
    fprintf(stdout,"Size=%u\n",ulCnt); 
}


void ProcessRecs()
{
   fprintf(stdout,"ProcessRecs...\n");

    VALIDITY_LIST::iterator itr; int nChroIdx=0; stKey *pstKey = new stKey;
    char acbuf[81920], acOut[81920]; char *pChr=NULL, *pcRsid=NULL;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; bool bFound;

    while (!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %uM\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'){
           fprintf(g_pfOut,"%s\n",acbuf); continue;}

        if (acbuf[0] == 'C'){fprintf(g_pfOut,"%s\tSNV_Validity\n",acbuf); continue;}

        strcpy(acOut,acbuf);

        pChr = strtok(acbuf,"\t"); //Chromosome  

        if (pChr[0]=='X') nChroIdx = 22;
        else if (pChr[0]=='Y') nChroIdx = 23; 
        else nChroIdx = atoi(pChr)-1;        

        pChr = strtok(NULL,"\t");  //GiNumber
        pChr = strtok(NULL,"\t");  //Offset
        pChr = strtok(NULL,"\t");  //Deleted_Base
        pChr = strtok(NULL,"\t");  //Fwd_DEL_Reads
        pChr = strtok(NULL,"\t");  //Rvs_DEL_Reads
        pChr = strtok(NULL,"\t");  //Total_DEL_Reads
        pChr = strtok(NULL,"\t");  //Total_Read_Density
        pChr = strtok(NULL,"\t");  //PEnd_Count
        pChr = strtok(NULL,"\t");  //Avg_QryPos
        pChr = strtok(NULL,"\t");  //Avg_QScore
        pChr = strtok(NULL,"\t");  //dbSNP

        fprintf(g_pfOut,"%s\t",acOut); bFound = false;
           
        if (pChr[0] != '-')
        {
            pcRsid = new char[strlen(pChr)+1]; strcpy(pcRsid, pChr);
            pChr = strtok(pcRsid,",");

            while(pChr)
            {
                if (pstKey->pcRsid) delete[] pstKey->pcRsid;

                pstKey->pcRsid = new char[strlen(pChr)+1];
                strcpy(pstKey->pcRsid,pChr);
 
                itr = g_Validity_List[nChroIdx].find(pstKey);

                if (itr != g_Validity_List[nChroIdx].end())
                {
                   fprintf(g_pfOut,"%s\n",(*itr->second).pcValidity);
                   bFound = true;
                   break;
                }
                pChr = strtok(NULL,",");
            }
        }

       if (!bFound){fprintf(g_pfOut,"-\n");}
    }
     
}


void CloseFiles()
{
    if (g_pfValidity) fclose(g_pfValidity);
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char **argv)
{
    if (argc < 4)
    {
       fprintf(stdout,"Invalid Parameters...\n");
       fprintf(stdout,"Usage:-\n"); 
       fprintf(stdout,"\t<SNV Validity filename> <Input filename> <Output filename> ");
       exit(9);   
    }

    g_pfValidity = fopen(argv[1],"r");

    if (!g_pfValidity){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfIn = fopen(argv[2],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }
 
    GenerateValidityList();ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); exit(9);
  
}
