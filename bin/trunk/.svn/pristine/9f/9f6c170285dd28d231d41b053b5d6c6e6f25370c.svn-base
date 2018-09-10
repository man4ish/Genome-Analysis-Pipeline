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
    char orient,*allele,*rsid;

    Value(const char* _rsid, const char _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1]; strcpy(rsid, _rsid);
        orient = _orient;
        allele = new char [strlen(_allele) + 1]; strcpy(allele, _allele);
    }

    ~Value()
    {
        if (rsid) delete [] rsid; if (allele) delete [] allele;
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->chronum != _rhs->chronum) return ( _lhs->chronum < _rhs->chronum);
        return  ( _lhs->offset < _rhs->offset );
    }
};


typedef multimap< Key*,Value*,ltstr> dbmap;

dbmap g_dbSNPMap; dbmap::iterator g_dbSNPit; pair<dbmap::iterator,dbmap::iterator> g_SNPret;
bool g_bFound,g_bdbSNP;

/*
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
*/

void  ReverseComplement(char *seq, char *pcRevSeq)
{
            
     int j=0;
     for (int i = strlen(seq)-1; i >= 0; i--)
     {
         switch (seq[i])
         {
                 case 'A': pcRevSeq[j++] = 'T'; break;
                 case 'C': pcRevSeq[j++] = 'G'; break;
                 case 'G': pcRevSeq[j++] = 'C'; break;
                 case 'T': pcRevSeq[j++] = 'A';// break;
        }
   }
}


int g_lendb; char *g_ndb, *g_pch1;

vector <string>g_dbvec;

bool match(const char allele ,const char orient,const char* db)
{
     g_dbvec.clear();

     if(db)
     {  
        g_lendb = strlen(db); g_ndb = new char[g_lendb+1]; strcpy(g_ndb,db);        
        g_pch1 = strtok(g_ndb,"/");

        while(g_pch1 != NULL)
        {
              //if (orient=='-') g_dbvec.push_back(ReverseComplement(g_pch1));              
              if (orient=='-'){
                  char *pcRevSeq = new char[strlen(g_pch1)+1];
                  ReverseComplement(g_pch1,pcRevSeq); g_dbvec.push_back(pcRevSeq);
                  delete[] pcRevSeq;    
              }
              else g_dbvec.push_back(g_pch1);
              
              g_pch1 = strtok(NULL,"/");
        }
        delete[] g_ndb;
     }

    for (unsigned int ii =0; ii< g_dbvec.size();ii++) { 
         if (g_dbvec[ii][0]==allele) return 1; 
    }

    return 0;
}


inline bool MapAllele( Key* location, const char cVariant, FILE* pf)
{   
     g_bFound = false;  
     g_SNPret = g_dbSNPMap.equal_range(location);
     for(g_dbSNPit=g_SNPret.first; g_dbSNPit!=g_SNPret.second; ++g_dbSNPit)
     {
         if (match(cVariant,(g_dbSNPit->second)->orient,(g_dbSNPit->second)->allele))
         {
             if (!g_bFound){ g_bFound = true; fprintf(pf,"%s",(g_dbSNPit->second)->rsid); }
             else { fprintf(pf,",%s",(g_dbSNPit->second)->rsid); }                         
         }
     }
     return g_bFound; 
}


bool OutputSNPNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{
    g_bdbSNP = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ()) {
        if (MapAllele(snpkey,cVariant,pf)) g_bdbSNP = true;
        else  fprintf(pf,"-");
    }
    else fprintf(pf,"-"); 

    return g_bdbSNP;
}


inline bool MapAllele( Key* location, const char cVariant, char acRsid[])
{
     g_bFound = false; g_SNPret = g_dbSNPMap.equal_range(location);
       
     for (g_dbSNPit=g_SNPret.first; g_dbSNPit!=g_SNPret.second; ++g_dbSNPit)
     {    
         if (match(cVariant,(g_dbSNPit->second)->orient,(g_dbSNPit->second)->allele))
         {
             if (!g_bFound){ g_bFound = true; sprintf(acRsid,"%s",(g_dbSNPit->second)->rsid); } 
             else { strcat(acRsid,","); strcat(acRsid,(g_dbSNPit->second)->rsid); }
         }

         if (acRsid[0]==0) sprintf(acRsid,"%s",(g_dbSNPit->second)->rsid); 
         else { strcat(acRsid,","); strcat(acRsid,(g_dbSNPit->second)->rsid); }                  
     }
     return g_bFound;
}


bool OutputSNPNovel(int chrno, size_t offset, const char cVariant,
                    FILE* pf, char acRsid[])
{
    g_bdbSNP = false; Key * snpkey = new Key(chrno, offset);   
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ()) {
        if (MapAllele(snpkey,cVariant,acRsid)) g_bdbSNP = true; else sprintf(acRsid,"-");        
    }
    else sprintf(acRsid,"-");             


    return g_bdbSNP;
}


bool GenerateDBSNPMap( FILE *dbfile)
{
    if (!dbfile) return false;

    char dbline[1000], *pch, *lpc[5]; int i, chr_num;

    while (fgets(dbline, 1000, dbfile))
    {
        i=0; pch = strtok (dbline,"\t");

        while(pch != NULL){
	  lpc[i++]=pch; if (i>4) break; pch = strtok (NULL,"\t");
        }

        if (lpc[1][0] == 'X') chr_num = 23;
        else if (lpc[1][0] == 'Y') chr_num = 24;        
        else chr_num = atoi(lpc[1]);
        
        Key * key = new Key(chr_num, (size_t)atol(lpc[2]));  // Key(int  _chronum, size_t  _offset):chronum(_chronum),offset(_offset)
	Value * value = new Value(lpc[0],lpc[3][0], lpc[4]); // Value(const char* _rsid, const char * _orient, const char * _allele)
	g_dbSNPMap.insert(pair<Key *, Value *>(key, value));

        memset(lpc,0,sizeof(lpc));
    }

    return true;
}


void ClrDBSNPMap()
{
    g_dbSNPMap.clear();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool MBmatch(const char *pcVar ,const char orient,const char* db)
{
     g_dbvec.clear();

     if(db)
     {
        g_lendb = strlen(db); g_ndb = new char[g_lendb+1]; strcpy(g_ndb,db);
        g_pch1 = strtok(g_ndb,"/");

        while(g_pch1 != NULL)
        {
              //if (orient=='-') g_dbvec.push_back(ReverseComplement(g_pch1));              
              if (orient=='-'){
                  char *pcRevSeq = new char[strlen(g_pch1)+1];
                  ReverseComplement(g_pch1,pcRevSeq); g_dbvec.push_back(pcRevSeq);
                  delete[] pcRevSeq;
              }
              else g_dbvec.push_back(g_pch1);

              g_pch1 = strtok(NULL,"/");
        }
        delete[] g_ndb;
     }

    for (unsigned int ii =0; ii< g_dbvec.size();ii++) {
         //if (g_dbvec[ii][0]==allele) return 1;
         if (strcmp(g_dbvec[ii].c_str(),pcVar)==0){
             //fprintf(stdout,"g_dbvec[ii].c_str() = %s,pcVar=%s\n",g_dbvec[ii].c_str(),pcVar);
             return 1; }
    }

    return 0;
}



inline bool MapMBAllele( Key* location, const char *pcVar, char acRsid[])
{       
     g_bFound = false; acRsid[0]=0; g_SNPret = g_dbSNPMap.equal_range(location);

     for (g_dbSNPit=g_SNPret.first; g_dbSNPit!=g_SNPret.second; ++g_dbSNPit)
     {
/*            
         if (MBmatch(pcVar,(g_dbSNPit->second)->orient,(g_dbSNPit->second)->allele))
         {
             if (!g_bFound){ g_bFound = true; sprintf(acRsid,"%s",(g_dbSNPit->second)->rsid); } 
             else { strcat(acRsid,","); strcat(acRsid,(g_dbSNPit->second)->rsid); }
         }
*/       

         g_bFound = true;
         if (acRsid[0]==0) sprintf(acRsid,"%s",(g_dbSNPit->second)->rsid);
         else { strcat(acRsid,","); strcat(acRsid,(g_dbSNPit->second)->rsid); }
       
     }
     return g_bFound;
}


bool OutputMBSNPNovel(int chrno, size_t offset, const char *pcVariant,
                    FILE* pf, char acRsid[])
{
    g_bdbSNP = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ()) {
       if (MapMBAllele(snpkey,pcVariant,acRsid)) g_bdbSNP = true; else sprintf(acRsid,"-");

       g_bdbSNP=true;

    }
    else sprintf(acRsid,"-");

    return g_bdbSNP;
}

