#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
#include <string>
using namespace std;

struct Key
{
    int gi;
    size_t offset;
    string allele;
    Key(int _gi, size_t  _offset, string _allele)
    {
      gi = _gi ;
      offset = _offset;
      allele = _allele; 
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi)
         return ( _lhs->gi < _rhs->gi);
        else if(_lhs->offset != _rhs->offset)
           return  ( _lhs->offset < _rhs->offset );
          return (strcmp((_lhs->allele).c_str(),(_rhs->allele).c_str())< 0);
    }
};

typedef map<Key*,string,ltstr> m;
m m1,m2,m3;
int commcount =0;
int fcount1 =0;
int fcount2 =0;
int reccount2  =0;
clock_t start = clock();
char Title[200];
int RecordSample1 = 0; 
int RecordSample2 = 0;

void SetTitle(const char* title)
{
   strcpy(Title,title);
}

const char* GetTitle()
{
   return Title;
}

void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && line1[0]!= '#') 
     {
           //if(line1[0] == 'C' && line1[1] == 'h' && line1[2] == 'r')
           //SetTitle(line1);
           char recline1[40000];
           strcpy(recline1,line1);
           char* pch;
           vector<string> sdf;
           string s;
           pch = strtok (recline1,"\t");
           while(pch != NULL) 
           {
              sdf.push_back(pch);
              pch = strtok (NULL,"\t");
           }
           int chrnum;
           size_t offset;
            if(sdf[0] == "X")
           chrnum = 23;
            else if(sdf[0]== "Y")
           chrnum = 24;
           else
           chrnum = atoi(sdf[0].c_str());
           Key* key= new Key(abs(chrnum),atoi(sdf[1].c_str()),"A"); //atoi(sdf[2].c_str()),"A");
           m1[key] = line1;
     }
}

void GetSNPs()
{
     for( m::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
     {
            
          if(strstr((ii->second).c_str(),"SNPS"))
          {
             bool flag = false;
             for (int i = -3 ; i < 4 ; i ++)
             {
                if(m1.count(new Key(ii->first->gi,ii->first->offset+i,"A")) != 0)
                {
                     if( strstr((m1.find(new Key(ii->first->gi,ii->first->offset+i,"A"))->second).c_str(),"INSS")|| strstr((m1.find(new Key(ii->first->gi,ii->first->offset+i,"A"))->second).c_str(),"DELS"))
                    flag = true; break;
                }
             }
             if(flag == true)
              cout << ii->second << endl;
          }
     }
}


int main(int argc ,char* argv[]) 
{
   string line1,line2;
   ifstream snpfile1 (argv[1]);
   if(!snpfile1.is_open())
   {
      cout<<"Unable to open the snpfile1\n";
      exit(1);
   }
   while(snpfile1) 
   {
          char line1[40000];
          snpfile1.getline(line1,40000);
          if(snpfile1)
            readsnpfile1(line1);
   }   
   snpfile1.close();
   cout << "snpfile1 loaded\n";
   GetSNPs();  
   m1.clear();
   return 0;
}
