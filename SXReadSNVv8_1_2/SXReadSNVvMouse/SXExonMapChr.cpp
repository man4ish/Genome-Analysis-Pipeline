/* 
 * File:   SXExonMapChr.cpp
 * Author: jameswong
 *
 * Created on October 21, 2009, 5:05 PM
 */

#include "SXExonMapChr.h"

EXON_LIST Exon_chrMap;
EXON_LIST::iterator Exon_chrIt;

bool GenerateExonMap(FILE *exonFile)
{
    if (!exonFile) return false;

    char qString[maxlen] = "\0";
    char * fields[20];

    size_t fCnt = 0;
    for (unsigned i = 0; i < 20; i++)
        fields[i] = NULL;

    size_t start,stop; unsigned chrNum;
    char *p,*chrPtr;
    while (fgets(qString, maxlen, exonFile))
    {
        p = strchr(qString, '\n');
        if (p != NULL)
            *p = '\0';
        if (qString[0] == '#')
            continue;
        fCnt = Exon_rowParser(qString, ",", fields, 20);
        chrPtr = strstr(fields[0], "NC_0") + 4;
        chrNum = (unsigned)atoi(chrPtr);

        start = (size_t)atol(fields[2] + 10);
        stop = (size_t)atol(fields[3] + 9);
        Exon_chrIt = Exon_chrMap.find(chrNum);
        if (Exon_chrIt == Exon_chrMap.end())
            Exon_chrIt = Exon_chrMap.insert(std::pair<unsigned, Exon_Chromosome*>(chrNum, new Exon_Chromosome())).first;

        Exon_chrIt->second->insertExon(start, stop, fields[5], fields[6]);
        for (unsigned i = 0; i < fCnt; i++)
            fields[i] = NULL;
    }
    
    for(Exon_chrIt = Exon_chrMap.begin();Exon_chrIt != Exon_chrMap.end();Exon_chrIt++)
    {
        Exon_chrIt->second->sortExon();
    }
    return true;
}


bool OutputExonID(unsigned int chrNum, unsigned int offset, FILE* pf)
{
    char * exonIds[20];
    Exon_chrIt = Exon_chrMap.find(chrNum);
    unsigned exonCnt = Exon_chrIt->second->search(offset, exonIds, 20);

    if (exonCnt != 0)
    {
        for (unsigned i = 0; i < exonCnt; i++)
        {
            if (i == 0)//exonCnt - 1)
                fprintf(pf,"%s", exonIds[i]);
            else
                fprintf(pf,",%s", exonIds[i]);
        }
    }
    else{
       fprintf(pf,"-"); return false;
    }
    return true;
}

void ClrExonMap()
{
    if (Exon_chrMap.empty() == false)
    {
        for (Exon_chrIt = Exon_chrMap.begin(); Exon_chrIt != Exon_chrMap.end(); Exon_chrIt++)
            delete Exon_chrIt->second;
        Exon_chrMap.clear();
    }
}

