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

using namespace std;

struct Key
{
    int gi;
    size_t offset;
    Key(int  _gi, size_t  _offset):gi(_gi),offset(_offset) {}
};

struct Value
{
    char orient, *allele, *rsid;

    Value(const char* _rsid, const char _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1]; strcpy(rsid, _rsid); 
        orient = _orient;
        allele = new char [strlen(_allele) + 1]; strcpy(allele, _allele);
    }

    ~Value()
    {
        if (rsid != NULL) delete [] rsid; if (allele != NULL) delete [] allele;                    
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi) return ( _lhs->gi < _rhs->gi);
        return  ( _lhs->offset < _rhs->offset );
    }
};

map<int ,long> ChrtoGi;
typedef multimap< Key*,Value*,ltstr> dbmap;

dbmap g_dbInDelMap; dbmap::iterator g_dbInDelit; pair< dbmap ::iterator,dbmap::iterator> g_InDelret;
bool g_bFound_INDEL, g_bdbINDEL;

/*
const char* ReverseComplementINDEL(char* seq)
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


void  ReverseComplementINDEL(char *seq, char *pcRevSeq)
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


int g_lendb_INDEL; char *g_ndb_INDEL, *g_pch1_INDEL;
vector <string>g_dbvec_INDEL;

bool matchINDEL(const char allele ,const char orient,const char* db)
{
     g_dbvec_INDEL.clear();

     if(db)
     {
        g_lendb_INDEL = strlen(db); g_ndb_INDEL = new char[g_lendb_INDEL+1]; strcpy(g_ndb_INDEL,db);        
        g_pch1_INDEL = strtok(g_ndb_INDEL,"/");

        while(g_pch1_INDEL != NULL)
        {
              //if (orient == '-') g_dbvec_INDEL.push_back(ReverseComplementINDEL(g_pch1_INDEL));
              if (orient=='-'){
                  char *pcRevSeq = new char[strlen(g_pch1_INDEL)+1];
                  ReverseComplementINDEL(g_pch1_INDEL,pcRevSeq); g_dbvec_INDEL.push_back(pcRevSeq);
                  delete[] pcRevSeq;          
              }
              else g_dbvec_INDEL.push_back(g_pch1_INDEL);
              
              g_pch1_INDEL = strtok(NULL,"/");
        }
        delete[] g_ndb_INDEL;
     }

     for (unsigned int ii =0; ii< g_dbvec_INDEL.size();ii++) {
          if (g_dbvec_INDEL[ii][0]==allele) return 1;      
     }
  
     return 0;
}


inline bool MapAllele( Key* location ,const char variant,  FILE* pf)
{
     g_bFound_INDEL = false;
     g_InDelret = g_dbInDelMap.equal_range(location);

     for (g_dbInDelit=g_InDelret.first; g_dbInDelit!=g_InDelret.second; ++g_dbInDelit)
     {    
         if (matchINDEL(variant,(g_dbInDelit->second)->orient,(g_dbInDelit->second)->allele))
         {
             if (!g_bFound_INDEL){ g_bFound_INDEL = true; fprintf(pf,"%s",(g_dbInDelit->second)->rsid); }
             else { fprintf(pf,",%s",(g_dbInDelit->second)->rsid); }
         }
     }
     return g_bFound_INDEL;
}


bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{    
    g_bdbINDEL = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbInDelMap.find (snpkey);

    if (elefound != g_dbInDelMap.end ())
    {
        if (MapAllele(snpkey,cVariant,pf)) g_bdbINDEL = true;        
        else fprintf(pf,"-");         
    }
    else fprintf(pf,"-");
    
    return g_bdbINDEL;
}


inline bool MapAllele( Key* location ,const char &variant, char acRsid[])
{
     g_bFound_INDEL = false; g_InDelret = g_dbInDelMap.equal_range(location);

     for (g_dbInDelit=g_InDelret.first; g_dbInDelit!=g_InDelret.second; ++g_dbInDelit)
     {
         if (matchINDEL(variant,(g_dbInDelit->second)->orient,(g_dbInDelit->second)->allele))
         {
             if (!g_bFound_INDEL){ g_bFound_INDEL = true; sprintf(acRsid,"%s",(g_dbInDelit->second)->rsid); }
             else { strcat(acRsid,","); strcat(acRsid,(g_dbInDelit->second)->rsid); }
         }
     }
     return g_bFound_INDEL;
}


bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf, char acRsid[])
{
    g_bdbINDEL = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbInDelMap.find (snpkey);

    if (elefound != g_dbInDelMap.end ()) {
        if (MapAllele(snpkey,cVariant,acRsid)) g_bdbINDEL = true;        
        else sprintf(acRsid,"-");        
    }
    else sprintf(acRsid,"-");    

    return g_bdbINDEL;
}


bool GenerateDBINDELMap(FILE *dbfile)
{
    if (!dbfile) return false;

    char dbline[1000], *pch, *lpc[5]; int i; int chr_num;

    while (fgets(dbline, 1000, dbfile))
    {
        i=0; pch = strtok (dbline,"\t");

        while(pch != NULL){
          lpc[i++]=pch; if (i>4) break; pch = strtok (NULL,"\t");
        }

        if (lpc[1][0] == 'X') chr_num = 23;        
        else if (lpc[1][0] == 'Y') chr_num = 24;        
        else chr_num = atoi(lpc[1]);
        
        Key * key = new Key(chr_num, (size_t)atol(lpc[2]));    // Key(int  _chronum, size_t  _offset):chronum(_chronum),offset(_offset)
        Value * value = new Value(lpc[0],lpc[3][0], lpc[4]);      // Value(const char* _rsid, const char * _orient, const char * _allele)
        g_dbInDelMap.insert(pair<Key *, Value *>(key, value));

        memset(lpc,0,sizeof(lpc));
    }

    return true;
}


void ClrDBINDELMap()
{
    g_dbInDelMap.clear();
}
