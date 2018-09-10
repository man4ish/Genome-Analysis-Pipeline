#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <iomanip>
#include <string.h>
#include "SXDBSNPgChecker_OdbSnpfmt.h"

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

dbmap g_dbSNPMap; dbmap::iterator g_dbSNPitr; pair<dbmap::iterator,dbmap::iterator> g_SNPretO;
bool g_bFoundO,g_bdbSNPO;

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


int g_lendbO; char *g_ndbO, *g_pchr1;

vector <string>g_dbvecO;

bool match(const char allele ,const char orient,const char* db)
{
     g_dbvecO.clear();

     if(db)
     {  
        g_lendbO = strlen(db); g_ndbO = new char[g_lendbO+1]; strcpy(g_ndbO,db);        
        g_pchr1 = strtok(g_ndbO,"/");

        while(g_pchr1 != NULL)
        {
              if (orient=='-'){
                  char *pcRevSeq = new char[strlen(g_pchr1)+1];
                  ReverseComplement(g_pchr1,pcRevSeq); g_dbvecO.push_back(pcRevSeq);
                  delete[] pcRevSeq;    
              }
              else g_dbvecO.push_back(g_pchr1);
              
              g_pchr1 = strtok(NULL,"/");
        }
        delete[] g_ndbO;
     }

    for (unsigned int ii =0; ii< g_dbvecO.size();ii++) { 
         if (g_dbvecO[ii][0]==allele) return 1; 
    }

    return 0;
}


inline bool MapAllele( Key* location, const char cVariant, FILE* pf)
{   
     g_bFoundO = false;  
     g_SNPretO = g_dbSNPMap.equal_range(location);
     for(g_dbSNPitr=g_SNPretO.first; g_dbSNPitr!=g_SNPretO.second; ++g_dbSNPitr)
     {
         if (match(cVariant,(g_dbSNPitr->second)->orient,(g_dbSNPitr->second)->allele))
         {
             if (!g_bFoundO){ g_bFoundO = true; fprintf(pf,"%s",(g_dbSNPitr->second)->rsid); }
             else { fprintf(pf,",%s",(g_dbSNPitr->second)->rsid); }                         
         }
     }
     return g_bFoundO; 
}


bool OutputSNPNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{
    g_bdbSNPO = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);


    if (elefound != g_dbSNPMap.end ()) {
        if (MapAllele(snpkey,cVariant,pf)) g_bdbSNPO = true;
        else  fprintf(pf,"-");
    }
    else fprintf(pf,"-"); 

    return g_bdbSNPO;
}


inline bool MapAllele( Key* location, const char cVariant, char acRsid[])
{
     g_bFoundO = false; g_SNPretO = g_dbSNPMap.equal_range(location);
       
     for (g_dbSNPitr=g_SNPretO.first; g_dbSNPitr!=g_SNPretO.second; ++g_dbSNPitr)
     {
         if (match(cVariant,(g_dbSNPitr->second)->orient,(g_dbSNPitr->second)->allele))
         {
             if (!g_bFoundO){ g_bFoundO = true; sprintf(acRsid,"%s",(g_dbSNPitr->second)->rsid); } 
             else { strcat(acRsid,","); strcat(acRsid,(g_dbSNPitr->second)->rsid); }
         }
     }
     return g_bFoundO;
}


bool OutputSNPNovel(int chrno, size_t offset, const char cVariant,
                    FILE* pf, char acRsid[])
{
    g_bdbSNPO = false; Key * snpkey = new Key(chrno, offset);   
    dbmap::const_iterator elefound = g_dbSNPMap.find (snpkey);

    if (elefound != g_dbSNPMap.end ()) {
        if (MapAllele(snpkey,cVariant,acRsid)) g_bdbSNPO = true; else sprintf(acRsid,"-");        
    }
    else sprintf(acRsid,"-");             

    return g_bdbSNPO;
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

