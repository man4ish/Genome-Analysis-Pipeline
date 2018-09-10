#include <string.h>
#include <iostream>
#include <map>
#include <utility>
#include <fstream>
#include <set>
#include <vector>
using namespace std;
multimap<int, string> m;
set<int> s;
vector<string> split(const char* record)
{
   char line[100000];
   strcpy(line,record);
   char* pch;
   vector<string> vec;
   pch = strtok(line,"\t");
   while(pch != NULL)
   { 
      vec.push_back(pch);
      pch = strtok(NULL,"\t");
   } 
   return vec;
}
int main(int argc, char* argv[])
{
  ifstream in(argv[1]);
  char line[100000];
  while(in)
  {
     in.getline(line,100000);
     if(in)
     {
        int chrnum ;
        string chr = split(line)[0];
        if(chr== "X")
        chrnum = 23;
        else  if(chr== "Y")
        chrnum = 24;
        else
        chrnum = atoi(chr.c_str());
        m.insert(pair<int, string>(chrnum,line));
        s.insert(chrnum);
     }
  }
  in.close();

   pair<multimap< int,string>::iterator, multimap< int,string>::iterator> ppp;
   set<int>::iterator it;
   for ( it=s.begin() ; it != s.end(); it++ )
   {
         ppp = m.equal_range(*it);
         multimap<int, string> m2;
         for (multimap<int,string>::iterator it2 = ppp.first;it2 != ppp.second;it2++)
         {
             string offset =   split(it2->second.c_str())[2];  // offset column index 
             m2.insert(pair<int, string>( atoi(offset.c_str()),(*it2).second ));
         }
         for (multimap<int,string>::iterator it = m2.begin();it != m2.end();it++)
         {
             cout<< (*it).second << endl;
         }
         m2.clear();
   }
   m.clear();
   s.clear(); 
}

