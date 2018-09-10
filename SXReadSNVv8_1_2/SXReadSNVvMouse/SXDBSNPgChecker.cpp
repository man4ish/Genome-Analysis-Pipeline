#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <iomanip>
#include "SXDBSNPgChecker.h"

using namespace std;

struct Key
{
    int chronum;
    size_t offset;
    Key(int  _chronum, size_t  _offset):chronum(_chronum),offset(_offset) {}
};

struct Value
{
    char * orient;
    char * allele ;
    char* rsid;
    Value(const char* _rsid, const char * _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1]; strcpy(rsid, _rsid);
        orient = new char [strlen(_orient) + 1]; strcpy(orient, _orient);
        allele = new char [strlen(_allele) + 1]; strcpy(allele, _allele);
    }

    ~Value()
    {
        if (rsid) delete [] rsid;
        if (orient) delete [] orient;
        if (allele) delete [] allele;
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->chronum != _rhs->chronum)
            return ( _lhs->chronum < _rhs->chronum);

        return  ( _lhs->offset < _rhs->offset );
    }
};


typedef multimap< Key*,Value*,ltstr> dbmap;
pair< multimap<Key*,Value*, ltstr> ::iterator,multimap<Key*,Value*, ltstr> ::iterator> g_SNPret;
dbmap g_dbSNPMap; dbmap::iterator g_dbSNPit;
int g_SNPnovcount =0,g_SNPknowncount =0,g_SNPcount =0;

const char* ReverseComplement(char* seq)
{
            string revseq;
            int l = strlen(seq);
            for( int i = l-1; i >= 0; i--)
            {
                 switch (seq[i])
                {
                         case 'A':
                         revseq += "T";
                         break;
                         case 'C':
                         revseq += "G";
                         break;
                         case 'G':
                         revseq += "C";
                         break;
                         case 'T':
                         revseq += "A";
                         break;
                }
            }
            return revseq.c_str();
}


int g_lendb; char *g_ndb, *g_pch1;

//char g_dbSNP_Alleele[1024];

vector <string>g_dbvec;

bool match(const char allele ,const char* orient,const char* db)
{
    g_dbvec.clear();

    //vector <string>g_dbvec;
     if(db)
     {  
        g_lendb = strlen(db);
        g_ndb = new char[g_lendb+1];
        strcpy(g_ndb,db);
        g_ndb[g_lendb] = '\0';

        g_pch1 = strtok(g_ndb,"/");
        while(g_pch1 != NULL)
        {
              if(strcmp(orient,"-")==0)
              {
                 g_dbvec.push_back(ReverseComplement(g_pch1));
              }
              else 
              {
                 g_dbvec.push_back(g_pch1);
              }
              g_pch1 = strtok(NULL,"/");
        }
        delete[] g_ndb;
     }

  for(unsigned int ii =0; ii< g_dbvec.size();ii++)
  {      
      if (g_dbvec[ii][0]==allele)
      {
          /*memset(g_dbSNP_Alleele,0,1024);
          for(unsigned int j =0; j< g_dbvec.size();j++)
          {

          }*/
         return 1;
      }
  }
  //g_dbvec.clear();
  return 0;
}


bool MapAllele( Key* location ,const char cVariant)
{     
     g_SNPret = g_dbSNPMap.equal_range(location);
     for(g_dbSNPit=g_SNPret.first; g_dbSNPit!=g_SNPret.second; ++g_dbSNPit)
     {
         if(match(cVariant,(g_dbSNPit->second)->orient,(g_dbSNPit->second)->allele)){
            return 1;
         }
     }
     return 0;
}



bool OutputSNPNovel2(int chrno, size_t offset,const char cVariant, FILE* pf)
{
    bool bdbSNP = false; g_SNPcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    short acACGT[4]={' ',' ',' ',' '};

    if (elefound != g_dbSNPMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
           //fprintf(pf,"Yes");
            fprintf(pf,"%s",(*elefound->second).rsid);

            for(unsigned int ii =0; ii< g_dbvec.size();ii++)
            {
                   //fprintf(pf,",%c",g_dbvec[ii][0]);
                if (g_dbvec[ii][0]=='A')
                    acACGT[0]='A';
                else if (g_dbvec[ii][0]=='C')
                    acACGT[1]='C';
                else if (g_dbvec[ii][0]=='G')
                    acACGT[2]='G';
                else if (g_dbvec[ii][0]=='T')
                    acACGT[3]='T';
            }

            for(unsigned int j =0; j< 4;j++)
            {
                fprintf(pf,"\t%c",acACGT[j]);
            }
            g_SNPknowncount++; bdbSNP = true;
        }
        else{
            fprintf(pf,"-\t \t \t \t");g_SNPnovcount++;
        }
    }
    else{
        fprintf(pf,"-\t \t \t \t"); g_SNPnovcount++;
    }
    return bdbSNP;
}


bool OutputSNPNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{
    bool bdbSNP = false; g_SNPcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
           //fprintf(pf,"Yes");
            fprintf(pf,"%s",(*elefound->second).rsid);
            g_SNPknowncount++; bdbSNP = true;
        }
        else{
            fprintf(pf,"-");g_SNPnovcount++;
        }
    }
    else{
        fprintf(pf,"-"); g_SNPnovcount++;
    }
    return bdbSNP;
}


bool OutputSNPNovel(int chrno, size_t offset, const char cVariant,
                    FILE* pf, char acRsid[])
{
    bool bdbSNP = false; g_SNPcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
            sprintf(acRsid,(*elefound->second).rsid);
            g_SNPknowncount++; bdbSNP = true;
        }
        else{
            sprintf(acRsid,"-");g_SNPnovcount++;
        }
    }
    else{
        sprintf(acRsid,"-"); g_SNPnovcount++;
    }
    return bdbSNP;
}


bool GenerateDBSNPMap( FILE *dbfile)
{
    if (!dbfile) return false;

    char dbline[1000]; char* pch;
    char *lpc[5];
    int i; int chr_num;

    while (fgets(dbline, 1000, dbfile))
    {
        i=0; pch = strtok (dbline,"\t");

        while(pch != NULL){

	  lpc[i++]=pch;
          if (i>4) break; pch = strtok (NULL,"\t");
        }

        if(strcmp(lpc[1],"X")==0){
           chr_num = 23;
        }
        else if(strcmp(lpc[1],"Y")==0){
           chr_num = 24;
        }
        else{
          chr_num = atoi(lpc[1]);
        }
        Key * key = new Key(chr_num, (size_t)atol(lpc[2]));    // Key(int  _chronum, size_t  _offset):chronum(_chronum),offset(_offset) {}
	Value * value = new Value(lpc[0],lpc[3], lpc[4]);      // Value(const char* _rsid, const char * _orient, const char * _allele)
	g_dbSNPMap.insert(pair<Key *, Value *>(key, value));

        memset(lpc,0,sizeof(lpc));
    }

    return true;
}


void GetSNPStatistics(const char* flname)
{
     ofstream out(flname);
     out<<"Total number of Snps : "<< g_SNPcount<<endl;
     out<<"Novel Snps           : "<<g_SNPnovcount<<"        "<<float (g_SNPnovcount*100)/float(g_SNPcount)<<"%"<<endl;
     out<<"Known Snps           : "<<g_SNPknowncount<<"      "<<float (g_SNPknowncount*100)/float (g_SNPcount)<<"%"<<endl;
     out.close();
}

void ClrDBSNPMap()
{
    g_dbSNPMap.clear();
}

