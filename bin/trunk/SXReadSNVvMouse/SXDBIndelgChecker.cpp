#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

struct Key
{
    int gi;
    size_t offset;
    Key(int  _gi, size_t  _offset):gi(_gi),offset(_offset) {}
};

struct Value
{
    char * orient;
    char * allele ;
    char* rsid;
    Value(const char* _rsid, const char * _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1];
        strcpy(rsid, _rsid); 
        orient = new char [strlen(_orient) + 1];
        strcpy(orient, _orient);
        allele = new char [strlen(_allele) + 1];
        strcpy(allele, _allele);
    }
    ~Value()
    {
        if (rsid != NULL)
        {
            delete [] rsid;
            rsid = NULL;
        }
 
        if (orient != NULL)
        {
            delete [] orient;
            orient = NULL;
        }
        if (allele != NULL)
        {
            delete [] allele;
            allele = NULL;
        }
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi)        
         return ( _lhs->gi < _rhs->gi);
           return  ( _lhs->offset < _rhs->offset );
    }
};

map<int ,long> ChrtoGi;
typedef multimap< Key*,Value*,ltstr> dbmap;
pair< multimap<Key*,Value*, ltstr> ::iterator,multimap<Key*,Value*, ltstr> ::iterator> g_InDelret;
dbmap g_dbInDelMap; dbmap::iterator g_dbInDelit;
int g_INDELnovcount =0,g_INDELknowncount =0,g_INDELcount =0;

/*
inline const char ReverseComplement(char seq)
{
    switch (seq)
    {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
    }

    return seq;
}


inline bool match(const char* allele ,const char* orient,const char* db)
{
      if(strcmp(orient,"-")== 0)
      {
         if(allele[0] == ReverseComplement(db[2]) )
         {
            return 1;
         } 
      }
      else 
      {
        if(allele[0] == db[2])
         {
            return 1;
         } 
      }   
  return 0;
}
*/


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



int g_lendb_INDEL; char *g_ndb_INDEL, *g_pch1_INDEL;
vector <string>g_dbvec_INDEL;

bool matchINDEL(const char allele ,const char* orient,const char* db)
{
    g_dbvec_INDEL.clear();

    //vector <string>g_dbvec;
     if(db)
     {
        g_lendb_INDEL = strlen(db);
        g_ndb_INDEL = new char[g_lendb_INDEL+1];
        strcpy(g_ndb_INDEL,db);
        g_ndb_INDEL[g_lendb_INDEL] = '\0';

        g_pch1_INDEL = strtok(g_ndb_INDEL,"/");
        while(g_pch1_INDEL != NULL)
        {
              //if(strcmp(&orient,"-")==0)
              if (strcmp(orient,"-")==0)
              {
                 g_dbvec_INDEL.push_back(ReverseComplementINDEL(g_pch1_INDEL));
              }
              else
              {
                 g_dbvec_INDEL.push_back(g_pch1_INDEL);
              }
              g_pch1_INDEL = strtok(NULL,"/");
        }
        delete[] g_ndb_INDEL;
     }

  for(unsigned int ii =0; ii< g_dbvec_INDEL.size();ii++)
  {
      if (g_dbvec_INDEL[ii][0]==allele)
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


inline bool MapAllele( Key* location ,const char &variant)
{
     
     g_InDelret = g_dbInDelMap.equal_range(location);
     for(g_dbInDelit=g_InDelret.first; g_dbInDelit!=g_InDelret.second; ++g_dbInDelit)
     {    
         if(matchINDEL(variant,(g_dbInDelit->second)->orient,(g_dbInDelit->second)->allele))
          {
            return 1;
          }
     }
     return 0;
}


bool GenerateDBINDELMap(FILE *dbfile)
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
	g_dbInDelMap.insert(pair<Key *, Value *>(key, value));

        memset(lpc,0,sizeof(lpc));
    }    

    return true;
}


bool OutputINDELNovel2(int chrno, size_t offset,const char cVariant, FILE* pf)
{    
    bool bdbSNP = false; g_INDELcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbInDelMap.find (snpkey);

    if (elefound != g_dbInDelMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
           //fprintf(pf,"Yes");
            fprintf(pf,"%s",(*elefound->second).rsid);
            g_INDELknowncount++; bdbSNP = true;
        }
        else{
            fprintf(pf,"-"); g_INDELnovcount++;
        }
    }
    else{
        fprintf(pf,"-"); g_INDELnovcount++;
    }
    return bdbSNP;
}


bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{    
    bool bdbSNP = false; g_INDELcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbInDelMap.find (snpkey);

    if (elefound != g_dbInDelMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
           //fprintf(pf,"Yes");
            fprintf(pf,"%s",(*elefound->second).rsid);
            g_INDELknowncount++; bdbSNP = true;
        }
        else{
            fprintf(pf,"-"); g_INDELnovcount++;
        }
    }
    else{
        fprintf(pf,"-"); g_INDELnovcount++;
    }
    return bdbSNP;
}


bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf, char acRsid[])
{
    bool bdbSNP = false; g_INDELcount++;
    Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = g_dbInDelMap.find (snpkey);

    if (elefound != g_dbInDelMap.end ())
    {
        if (MapAllele(snpkey,cVariant)==1){
            sprintf(acRsid,(*elefound->second).rsid);
            g_INDELknowncount++; bdbSNP = true;
        }
        else{
            sprintf(acRsid,"-"); g_INDELnovcount++;
        }
    }
    else{
        sprintf(acRsid,"-"); g_INDELnovcount++;
    }

    return bdbSNP;
}


void GetINDELStatistics(const char* flname)
{
     ofstream out(flname);
     out<<"Total number of Indels : "<< g_INDELcount<<endl;
     out<<"Novel Indels           : "<<g_INDELnovcount<<"        "<<float (g_INDELnovcount*100)/float(g_INDELcount)<<"%"<<endl;
     out<<"Known Indels           : "<<g_INDELknowncount<<"      "<<float (g_INDELknowncount*100)/float (g_INDELcount)<<"%"<<endl;
     out.close();
}


void ClrDBINDELMap()
{
    g_dbInDelMap.clear();
}
