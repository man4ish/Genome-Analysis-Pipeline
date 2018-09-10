/* 
 * File:   Avg_Confi_Score.cpp
 * Author: wuntyng
 *
 * Created on June 16, 2011, 3:06 PM
 */

#include <cstdlib>
#include <string.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include <map>

using namespace std;

typedef unsigned UINT;

inline size_t chomp(char *_string)
{
    size_t len = strlen(_string);
    while (_string[len - 1] == '\n')
        _string[--len] = '\0';
    while (_string[len - 1] == '\r')
        _string[--len] = '\0';
    return len;
}

inline size_t rowParser(char * _source, const char * _delimiter, char ** _fields, size_t _max)
{
    size_t index = 0;
    char * pch = strtok(_source, _delimiter);
    while((pch != NULL) && (index < _max))
    {
        _fields[index++] = pch;
        pch = strtok(NULL, _delimiter);
    }
    return index;
}

void manual()
{
    fprintf(stderr, "\n  Usage: Avg_Confi_Score <score file> <summary file> <output file>\n\n");
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        manual();
        exit(0);
    }
    FILE * inFile_score = NULL;
    FILE * inFile_summary = NULL;
    FILE * outFile = NULL;
    char lines[40960] = {0};
    char qString[2000] = {0};
    char *fields[100];
    int cnt = 0;
    int sr = 0;
    double cs = 0.000;
    std::map<int, double> sup_read_map;
    std::map<int, int> conf_map;
    std::map<int, double>::iterator sup_read_it;
    std::map<int, int>::iterator conf_it;

    inFile_score = fopen(argv[1], "r");
    if (inFile_score == NULL)
    {
        printf("Error: Unable to open the file %s\n", argv[1]);
        return false;
    }

    inFile_summary = fopen(argv[2],"r");
    if (inFile_summary == NULL)
    {
        printf("Error: Unable to open the file %s\n", argv[2]);
        exit(1);
    }

    while (fgets(lines, 40960, inFile_score))
    {
        chomp(lines);
        cnt = rowParser(lines, "\t", fields, 100);

        if(fields[0][0] == 'C'||fields[0][0]=='#')
            continue;
        
        sr = (int)atoi(fields[6]);
        cs = (double)atof(fields[cnt-1]);

        sup_read_it = sup_read_map.find(sr);
        if (sup_read_it == sup_read_map.end())
            sup_read_map.insert(std::pair<int, double>(sr, cs));
        else
            sup_read_it->second+=cs;

        conf_it = conf_map.find(sr);
        if (conf_it == conf_map.end())
            conf_map.insert(std::pair<int, int>(sr, 1));
        else
            conf_it->second++;
    }
    fclose(inFile_score);

//    Checking purpose
//    for (sup_read_it = sup_read_map.begin(); sup_read_it != sup_read_map.end(); sup_read_it++)
//    {
//        printf("%d\t%f\n", sup_read_it->first, sup_read_it->second);
//    }
//    for (conf_it = conf_map.begin(); conf_it != conf_map.end(); conf_it++)
//    {
//        printf("%d\t%d\n", conf_it->first, conf_it->second);
//    }

    outFile = fopen(argv[3],"w");
    if (outFile == NULL)
    {
        printf("Error: Unable to open the file %s\n", argv[3]);
        exit(1);
    }

    int fnt = 0; char *pChr=NULL;
    while (fgets(qString, 2000, inFile_summary))
    {
        chomp(qString);

        if(qString[0] != '=')
        {
            if (qString[0]!='S' && qString[0]!='>')
               fprintf(outFile, "%s\n", qString);
            else
            {
               pChr=strtok(qString,"\t"); fprintf(outFile, "%s", pChr);  

               for (int i=0; i <8; i++) {
                   pChr=strtok(NULL,"\t"); fprintf(outFile, "\t%s", pChr);
               }
              
               fprintf(outFile, "\t");
               if (qString[0] =='S') fprintf(outFile, "Avg_Confident_Score");  
               fprintf(outFile, "\t");  

               for (int i=0; i <9; i++) {
                   pChr=strtok(NULL,"\t"); fprintf(outFile, "\t%s", pChr);
               }

               fprintf(outFile,"\n");                 
            }
            continue;
        }

        fnt = rowParser(qString, "\t", fields, 100);

        //print the first half
        for(int i = 0; i < 9; i++)
            fprintf(outFile, "%s\t", fields[i]);

        char *seq_ptr = strstr(fields[0], "=");
        int sup_read = atoi(seq_ptr + 1);
        double avg = sup_read_map[sup_read]/conf_map[sup_read];
        if(avg > 0)
            fprintf(outFile, "%.4f\t\t", avg);
        else
            fprintf(outFile, "0.0\t\t");

        //print the rest
        for(int i = 9; i < fnt; i++)
            fprintf(outFile, "%s\t", fields[i]);
        fprintf(outFile, "\n");

    }
    fclose(inFile_summary);
    fclose(outFile);

    return 0;
}
