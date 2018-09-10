/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on May 10, 2010, 1:56 PM
 */

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string.h>
#include <ctype.h>
/*
 * 
 */

using namespace std;
map <int,int> m_List;

int main(int argc, char** argv) {

    
    if (argc < 3) //4)
    {
        fprintf(stdout,"Invalid parameters...\n");
        fprintf(stdout,"<Input filename> <Output filename>\n");
        //fprintf(stdout,"<Input filename> <Type-S,I,D> <Output filename>\n");   
        exit(9);
    }    
    
    FILE *pfIn = fopen(argv[1],"r");
    if(!pfIn){fprintf(stdout,"Failed to open %s ...\n",argv[1]);exit(9);}

    /*
    char cSNVType = toupper(argv[2][0]);
    if (cSNVType!= 'S' && cSNVType!= 'I' && cSNVType != 'D')
       {fprintf(stdout,"Invalid Type Specified %c ...\n",cSNVType);exit(9);}            

    FILE *pfOut = fopen(argv[3],"w");     
    if(!pfOut){fprintf(stdout,"Failed to open %s ...\n",argv[3]); exit(9);}
    */

    FILE *pfOut = fopen(argv[2],"w");
    if(!pfOut){fprintf(stdout,"Failed to open %s ...\n",argv[2]); exit(9);}


    //FILE *pfChk = fopen("ChkRD.lst","w");
 
    #if _TESTING
        unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;
    #endif

    
    char acbuf[1024]; char *pChr=NULL; char cChromosome;
    int nOffset, nSnpCnt, nRS;
    unsigned int unRD;
    
    while (!feof(pfIn))
    {
        if (!fgets(acbuf,sizeof(acbuf),pfIn)) break;
        
        pChr=strchr(acbuf,'\r'); if (pChr) *pChr=0;
        pChr=strchr(acbuf,'\n'); if (pChr) *pChr=0;
        
        if (acbuf[0] == '#'||acbuf[0] == 'C'||acbuf[0] == 0){         
           continue;
        }
        
        pChr = strtok(acbuf,"\t");   //chromosome

        if (!isalpha(pChr[0])) cChromosome=atoi(pChr);
        else if (pChr[0]=='X') cChromosome=23;
        else if (pChr[0]=='Y') cChromosome=24;
        else {fprintf(stdout,"Invalid chromomose number %d ...", atoi(pChr)); continue;}
        
        pChr = strtok(NULL,"\t");   //offset
        nOffset = atoi(pChr);
        
        pChr = strtok(NULL,"\t");   //Allele
        pChr = strtok(NULL,"\t");   //Total SNP_Density
        nSnpCnt = atoi(pChr);

        pChr = strtok(NULL,"\t");   //Total Read_Density
        unRD = atoi(pChr);  

        /* 
        if (cSNVType == 'I'){
            nDiff = unRD - nSnpCnt;
            if (nDiff > nSnpCnt) nSnpCnt=nDiff;
        } 
        */        
        
        nRS=(int)((double)(nSnpCnt*100)/(double)(unRD));

        if (nRS > 100 || nRS < 0) nRS = 101;
        
        if (m_List.count(nRS)==0) m_List[nRS] = 1;
        else m_List[nRS]++;
    }

//fclose(pfChk);

    for (map<int,int>::iterator itr = m_List.begin(); itr != m_List.end(); itr++)
        fprintf(pfOut,"%d\t%d\n",itr->first,itr->second);

    if (pfIn) fclose(pfIn); if (pfOut) fclose(pfOut);            
    
    return (EXIT_SUCCESS);
}

