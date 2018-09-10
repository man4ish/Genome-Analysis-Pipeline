/* 
 *
 * File:   main.cpp
 * Author: Manish
 *
 * Created on December 9, 2009, 4:23 PM
 */

#include "SSNpParser.h"
/*
 * 
 */

#define TOTAL_QRY_FIELDS 5
#define VARIANT_POS     45

FILE *g_pfSnpLst, *g_pfOut, *g_pfIn,*g_pfStats;
char g_acBuf[40960], *g_pChr;
const char *g_pcOutDir;
unsigned int g_unBufSize = sizeof(g_acBuf),g_unTotalExons;
SYNO_LIST g_SynoList; stStats g_stStats;


static void banner(char *argv[])
{
    printf("Synamatix SSNpParser Copyright 2009 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s", argv[0]);
    printf(" <Input file generated by SynaSearch> <SNV_SNP list filename> <output filename> <Sample_id> <output dir-stats. files>\n\n");
}


void ProcessSSOutFile()
{
    char * fields[4]; char acQrySeq[100],acSbjSeq[100]; 
    int nIdx=0, nQryStart=0, nQryStop=0, nSbjStart=0, nSbjStop=0;
    bool bGeneFound=true,bInRange=true, bCompleted=false, bNoHit=false;
    float fCodon; unsigned unSeqIdx;
    stSyno *pstSyno=NULL; stKey *pKey=NULL;

    while (!feof(g_pfIn))
    {
        if (!fgets(g_acBuf,g_unBufSize,g_pfIn)) break;
        g_pChr = strchr(g_acBuf, '\n'); if (g_pChr) {*g_pChr = '\0';}
        g_pChr = strchr(g_acBuf, '\r'); if (g_pChr) {*g_pChr = '\0';}
        if ((g_acBuf[0] == '\0') || (g_acBuf[0] == '#')) {continue;}
        
        if ((strncmp(g_acBuf, "Query=", 6) == 0))
        {
            if (bNoHit){      
                memcpy(pstSyno->acStatus,"NH",2); bCompleted = true;
            }
            else if (!bGeneFound){                
                memcpy(pstSyno->acStatus,"MM",2); bCompleted = true;
            }
            else if (!bInRange){                
                memcpy(pstSyno->acStatus,"OR",2); bCompleted = true;
            }

            if (bCompleted)
                g_SynoList[pKey]=pstSyno;
            
            bGeneFound=bInRange=bCompleted=bNoHit=false; nIdx=0;                       
            pKey = new stKey;  pstSyno = new stSyno;

            g_pChr = strtok(g_acBuf+7,"|");

            while (g_pChr){
                if (nIdx==1)
                    pKey->unGiNum = atol(g_pChr);
                else if (nIdx==2)
                    pKey->unOffset = atol(g_pChr);
                else if (nIdx==3)
                {
                    pKey->cRefBase = g_pChr[0];pKey->cVarBase = g_pChr[1];
                }
                else if (nIdx==4){
                    pstSyno->pcGeneIDs = new char[strlen(g_pChr)+1];
                    strcpy(pstSyno->pcGeneIDs,g_pChr);
                    strupr(pstSyno->pcGeneIDs);
                }
                if (++nIdx == TOTAL_QRY_FIELDS) break;
                g_pChr=strtok(NULL,"|");
            }            
        }
        else if (!bCompleted && g_acBuf[0] == '>')
        {
            g_pChr = strstr(g_acBuf,"OS=Homo sapiens"); if (!g_pChr) continue;
            g_pChr = strstr(g_acBuf,"GN="); if (!g_pChr) continue;
            strupr(g_pChr); 

            if (strstr(g_pChr,pstSyno->pcGeneIDs)) bGeneFound = true;
        }

        else if (!bCompleted && bGeneFound && strncmp(g_acBuf, "Query:", 6) == 0)
        {
            rowParser(g_acBuf, " ", fields, 4);
            nQryStart = (unsigned short)atoi(fields[1]);
            nQryStop = (unsigned short)atoi(fields[3]);
            strcpy(acQrySeq, fields[2]);

            bInRange = inRange(nQryStart, nQryStop, VARIANT_POS);
            if (!bInRange) bGeneFound=false;
        }
        else if (!bCompleted && bGeneFound && bInRange && strncmp(g_acBuf, "Sbjct:", 6) == 0)
        {
            rowParser(g_acBuf, " ", fields, 4); strcpy(acSbjSeq, fields[2]);
            nSbjStart = (size_t)atoi(fields[1]); nSbjStop = (size_t)atoi(fields[3]);

            if (nQryStart < nQryStop){
                fCodon = (VARIANT_POS - nQryStart + 3)/3.0;
                unSeqIdx = (unsigned)fCodon -1;
            }
            else{
                fCodon = (VARIANT_POS - nQryStop + 3)/3.0;
                unSeqIdx = strlen(acQrySeq) - (unsigned)fCodon;
            }            

            if (acSbjSeq[unSeqIdx] != acQrySeq[unSeqIdx]){
                memcpy(pstSyno->acStatus,"NS",2);         
                pstSyno->acProtVar[0] = acSbjSeq[unSeqIdx];
                pstSyno->acProtVar[1] = acQrySeq[unSeqIdx];
                
            }
            else                
                memcpy(pstSyno->acStatus,"S ",2);

            bCompleted = true;
        }
        else if (strncmp(g_acBuf, "No significant", 14) == 0)
            bCompleted = bNoHit = true;        
    }

    if (!bGeneFound)        
        memcpy(pstSyno->acStatus,"MM",2);    
    else if (!bInRange)        
        memcpy(pstSyno->acStatus,"OR",2);
        
    g_SynoList[pKey]=pstSyno;

/*    fprintf(stdout,"total =%u\n",g_SynoList.size());

    FILE *pfO2 = fopen("O2.tst","w");
    unsigned long ulcnt=0;

    for (SYNO_LIST::iterator itr=g_SynoList.begin(); itr!=g_SynoList.end(); itr++)
    {
        ulcnt++;
        if (ulcnt > 500 ) break;
        fprintf(pfO2,"%u->%u|%u|%c%c|%s|%c%c|\n",ulcnt,
                (*itr->first).unGiNum,
                (*itr->first).unOffset,
                (*itr->first).cRefBase,(*itr->first).cVarBase,
                (*itr->second).pcGeneIDs,(*itr->second).acStatus[0],(*itr->second).acStatus[1]);
    }
    fclose(pfO2);
*/
}


void PrintNewSNPList()
{
    unsigned int nCnt=0, nCol_Exon=0; SYNO_LIST::iterator itr; stKey SynoKey;

    while (!feof(g_pfSnpLst))
    {
        if (!fgets(g_acBuf, g_unBufSize,g_pfSnpLst)) break;
        
        g_pChr = strchr(g_acBuf,'\n'); if (g_pChr) {*g_pChr='\0';}           
        g_pChr = strchr(g_acBuf,'\r'); if (g_pChr) {*g_pChr='\0';}  

        if (g_acBuf[0] == '\n'){
           fprintf(g_pfOut,"\n"); continue;
        }

        if (g_acBuf[0] == '#') { 
           fprintf(g_pfOut,"%s\n",g_acBuf); continue;
        } 

        if (g_acBuf[0] == 'C'){
           fprintf(g_pfOut,"%s\n",g_acBuf); g_pChr = strtok(g_acBuf,"\t");

           while(g_pChr){
                 nCnt++; 
                 if (strcmp(g_pChr,"Exon")==0){
                     nCol_Exon=nCnt; break;
                 }
                 g_pChr = strtok(NULL,"\t");
           }

           if (nCol_Exon != 0) continue;
           else {fprintf(stdout, "Failed to local Exon fields !!!"); return;}  
        }

        fprintf(g_pfOut,"%s",g_acBuf);

        g_pChr = strchr(g_acBuf,'\t'); nCnt=1;

        while(g_pChr){
	    if(++nCnt==nCol_Exon) break; //syhwah
            g_pChr = strchr(g_pChr+1,'\t');
        }

        if (!g_pChr) continue; 

        if (g_pChr[1]=='-'){
            fprintf(g_pfOut,"\t-\t-\t-\t-\n"); continue;
        }
        else{
            g_unTotalExons++;

            nCnt=1; g_pChr = strtok(g_acBuf,"\t");
            while (g_pChr)
            {
                if (nCnt==2)
                    SynoKey.unGiNum= atol(g_pChr);
                else if (nCnt==3)
                    SynoKey.unOffset = atol(g_pChr);
                else if (nCnt==4){
                    SynoKey.cRefBase = g_pChr[0];SynoKey.cVarBase = g_pChr[2];
                }
                if (++nCnt > 4) break;  g_pChr = strtok(NULL,"\t");
            }

            itr = g_SynoList.find(&SynoKey);

            if (itr!=g_SynoList.end()){

#if _TESTING                 
            if (ulCnt > 500) break;
#endif                
                g_stStats.unTotalRec++; 

                if ((*itr->second).acStatus[0]=='N'&&(*itr->second).acStatus[1]=='S')
                {

                    fprintf(g_pfOut,"\tNS\t%c>%c",
                            (*itr->second).acProtVar[0],(*itr->second).acProtVar[1]);

                    if ((*itr->second).acProtVar[1]!='*')
                        fprintf(g_pfOut,"\tYes\tNo");
                    else
                        fprintf(g_pfOut,"\tNo\tYes");

                    g_stStats.unNS++;
                }
                else
                {
                    if ((*itr->second).acStatus[0]=='S'){
                        fprintf(g_pfOut,"\tS");g_stStats.unS++;
                    }
                    else{
                        fprintf(g_pfOut,"\t-");
                        if ((*itr->second).acStatus[0]=='N'&&(*itr->second).acStatus[1]=='H')
                            g_stStats.unNH++;
                        else if ((*itr->second).acStatus[0]=='M'&&(*itr->second).acStatus[1]=='M')
                            g_stStats.unMM++;
                        else if ((*itr->second).acStatus[0]=='O'&&(*itr->second).acStatus[1]=='R')
                            g_stStats.unOR++;
                    }
                    fprintf(g_pfOut,"\t-\t-\t-");
                }

                fprintf(g_pfOut,"\n");
            }
            else{
                fprintf(g_pfOut,"\t-\t-\t-\t-\n");  
            }
        }
    }
}


void PrintStatsRpt(const char* pcSampleID)
{
    fprintf(g_pfStats,"#REPORT NAME\tSynonymous/Non-Synonymous Statistic\n");
    fprintf(g_pfStats,"#PROJECT NAME\n");
    fprintf(g_pfStats,"#SAMPLE ID\t%s\n",pcSampleID);
    fprintf(g_pfStats,"#LANE NO\tALL\n");
    fprintf(g_pfStats,"#GENERATED AT\t"); PrintRptDateTime(g_pfStats);
    fprintf(g_pfStats,"#PROGRAM & BUILD\tSSNpParser Rev 1.0.0\n");
    fprintf(g_pfStats,"#REMARKS\t");PrintRemarksParam(g_pfStats);
    fprintf(g_pfStats,"#FILTER\t>=1\n\n\n");

    fprintf(g_pfStats,"\tTotal\tPercentage(%%)\n");
    fprintf(g_pfStats,"Total SNPs in Exon\t%u\n",g_stStats.unTotalRec);
    fprintf(g_pfStats,"Synonymous(S)\t%u\t%.2f\n",g_stStats.unS,Percentage(g_stStats.unS,g_stStats.unTotalRec));
    fprintf(g_pfStats,"Non-Synonymous(NS)\t%u\t%.2f\n",g_stStats.unNS,Percentage(g_stStats.unNS,g_stStats.unTotalRec));
    fprintf(g_pfStats,"Missed Matched\t%u\t%.2f\n",g_stStats.unMM,Percentage(g_stStats.unMM,g_stStats.unTotalRec));
    fprintf(g_pfStats,"No Hit\t%u\t%.2f\n",g_stStats.unNH,Percentage(g_stStats.unNH,g_stStats.unTotalRec));
    fprintf(g_pfStats,"Out of Range\t%u\t%.2f\n",g_stStats.unOR,Percentage(g_stStats.unOR,g_stStats.unTotalRec));

    fprintf(g_pfStats,"\n\nTotal Exons=\t%u",g_unTotalExons); 
     
}


int main(int argc, char** argv) {
    if (argc <6) {banner(argv); exit(9);}

    struct tm *ptm; time_t ttime;
    (void) time(&ttime); ptm = localtime(&ttime);

    char acRptName[2048];

    g_pfIn = fopen(argv[1], "r");
    if (!g_pfIn) {fprintf(stderr,"Failed to open %s ....\n",argv[1]); goto Exit;}

    g_pfSnpLst = fopen(argv[2], "r");
    if (!g_pfSnpLst) {fprintf(stderr,"Failed to open %s ....\n",argv[2]); goto Exit;}

    g_pfOut = fopen(argv[3], "w");
    if (!g_pfOut) {fprintf(stderr,"Failed to open %s ....\n",argv[3]); goto Exit;}
     
    g_pcOutDir = argv[5]; 

    sprintf(acRptName,"%s/%s_syno_stat_%02d%02d%02d.rpt",g_pcOutDir,
                      argv[4],ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

    g_pfStats = fopen(acRptName,"w");
    if (!g_pfStats) {fprintf(stderr,"Failed to open %s ....\n",acRptName); goto Exit;}

    ProcessSSOutFile(); 


    if (g_SynoList.size()){
        PrintNewSNPList(); PrintStatsRpt(argv[4]); 
    }
 Exit:    
    if (g_pfIn) fclose(g_pfIn); if (g_pfSnpLst) fclose(g_pfSnpLst); 
    if (g_pfOut) fclose(g_pfOut); if (g_pfStats) fclose(g_pfStats);

    SYNO_LIST::iterator itr = g_SynoList.begin();
    
    while (itr!= g_SynoList.end())
    {
       delete itr->first; delete itr->second;
       itr++;
    }
    g_SynoList.clear();
    return (EXIT_SUCCESS);
}

    


