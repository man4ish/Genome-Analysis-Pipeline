/* 
 * File:   SXCNVMap.cpp
 * Author: Manish
 * 
 * Created on February 2, 2010, 5:13 PM
 */

#include "SXCNVMap.h"

std::map<unsigned, CNV_Chromosome*> CNV_chrMap;
std::map<unsigned, CNV_Chromosome*>::iterator CNV_chrIt;

bool GenerateCNVMap(FILE *pf)
{
    if (!pf) return false;    

    char qString[1024] = "\0";
    char * fields[20];
    int maxlen = sizeof(qString);
       
    size_t start,stop; unsigned chrNum;

    while (fgets(qString, maxlen, pf))
    {
        if (qString[0] == '#')
            continue;
        CNV_rowParser(qString, "\t", fields, 20);

        if (fields[2][3] == 'X')
            chrNum = 23;
        else if (fields[2][3] == 'Y')
            chrNum = 24;
        else
            chrNum = (unsigned)atoi(&fields[2][3]);

        start = (size_t)atol(fields[3]);
        stop = (size_t)atol(fields[4]);
        CNV_chrIt = CNV_chrMap.find(chrNum);

        if (CNV_chrIt == CNV_chrMap.end())
            CNV_chrIt = CNV_chrMap.insert(std::pair<unsigned, CNV_Chromosome*>(chrNum, new CNV_Chromosome())).first;

        CNV_chrIt->second->insertCNV(start, stop);
    }

/*     
    char acbuf[15];
    for (int i=1;i<25;i++)
    {
        CNV_chrIt = CNV_chrMap.find(i); sprintf(acbuf,"CNV_%02d.lst",i);
        //fprintf(stdout,"%size=u\n",CNV_chrIt.);

        if (CNV_chrIt != CNV_chrMap.end())
            CNV_chrIt->second->printAll(acbuf);
    }
*/    

    return true;
}


void OutputCNVSearch(unsigned int chrNum, unsigned int offset, FILE* pf)
{
    CNV_chrIt = CNV_chrMap.find(chrNum);

    if (CNV_chrIt==CNV_chrMap.end()){
        fprintf(pf,"N"); return;
    }

    if (CNV_chrIt->second->search(offset))
        fprintf(pf,"Y");
    else
        fprintf(pf,"N");
}

void ClrCNVMap()
{
    if (CNV_chrMap.empty() == false)
    {
        for (CNV_chrIt = CNV_chrMap.begin(); CNV_chrIt != CNV_chrMap.end(); CNV_chrIt++)
            delete CNV_chrIt->second;
        CNV_chrMap.clear();
    }    
}
