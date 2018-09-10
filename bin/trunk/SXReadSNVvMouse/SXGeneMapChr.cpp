/* 
 * File:   SXGeneMapChr.cpp
 * Author: jameswong
 * Modified/Enhanced by: syhwah
 *
 * Created on October 21, 2009, 4:52 PM
 */

#include "SXGeneMapChr.h"
#define maxlen 1024


std::map<unsigned, Chromosome*> chrMap;
std::map<unsigned, Chromosome*>::iterator chrIt;

bool GenerateGeneMap(FILE *geneFile)
{
    if (!geneFile) return false;

    char qString[maxlen] = "\0";
    char * fields[20]; char *pChr;

    size_t fCnt,start,stop; unsigned chrNum;

    while (fgets(qString, maxlen, geneFile))
    {
        if (qString[0] == '#') continue;

        pChr = strchr(qString,'\n'); if (pChr) *pChr=0;
        pChr = strchr(qString,'\n'); if (pChr) *pChr=0;  
          
        fCnt = rowParser(qString, "\t", fields, 20);
        chrNum = 0;
        if (fields[2][0] == 'X')
            chrNum = 23;
        else if (fields[2][0] == 'Y')
            chrNum = 24;
        else
            chrNum = (unsigned)atoi(fields[2]);

        start = (size_t)atol(fields[3]);
        stop = (size_t)atol(fields[4]);
        chrIt = chrMap.find(chrNum);
        if (chrIt == chrMap.end())
            chrIt = chrMap.insert(std::pair<unsigned, Chromosome*>(chrNum, new Chromosome())).first;
        chrIt->second->insertGene(start, stop, fields[1],fields[6]);
    }    
    return true;
}

char * g_geneIds[20];char * g_geneDescs[20];unsigned g_geneCnt; bool g_bGene;
bool OutputGeneIDnDesc(unsigned int chrNum, unsigned int offset, FILE* pf)
{
    //fprintf(pf,"-\t-"); return true;

    for (int i=0;i<20;i++)
    {
        if (!g_geneIds[i]) break;
        g_geneIds[i]=NULL; g_geneDescs[i]=NULL;
    }

    /*for (int i=0;i<20;i++)
    {
        if (!g_geneDescs[i]) break;
        g_geneDescs[i]=NULL;
    }*/

    g_bGene = true; chrIt = chrMap.find(chrNum);
    g_geneCnt = chrIt->second->search(offset, g_geneIds,g_geneDescs, 20);

    if (g_geneCnt != 0)
    {
       for (unsigned i = 0; i < g_geneCnt; i++)
       {
            if (i == 0)//g_geneCnt - 1)
                fprintf(pf,"%s", g_geneIds[i]);
            else
                fprintf(pf,";%s", g_geneIds[i]);
       }
       fprintf(pf,"\t");

       bool bNA =false;

       for (unsigned i = 0; i < g_geneCnt; i++)
       {
            if (strcmp(g_geneDescs[i],"NA")==0){
               if (!bNA) bNA = true;
               else continue;
            }
            
            if (i == 0){
                if (strlen(g_geneDescs[i])<2) fprintf(pf,"-");  
                else fprintf(pf,"%s", g_geneDescs[i]);
            }
            else{
                if (strlen(g_geneDescs[i])<2) fprintf(pf,"-"); 
                else fprintf(pf,";%s", g_geneDescs[i]);
             }
       }
    }
    else
    {
       fprintf(pf,"-\t-"); g_bGene = false;

       if (!g_geneIds[0])
           g_geneIds[0] = new char[2];

       strcpy(g_geneIds[0],"-");
    }


    return g_bGene;
}


void ClrGeneMap()
{
    if (chrMap.empty() == false)
    {
        for (chrIt = chrMap.begin(); chrIt != chrMap.end(); chrIt++)
            delete chrIt->second;
        chrMap.clear();
    }    
}

