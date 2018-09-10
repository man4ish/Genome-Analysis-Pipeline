#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
map <int,int> m; 
void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chromosome"))
           return;
        else
        {
           char recline1[10000];
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
           int snp =  atoi(sdf[5].c_str());
           int diff =  atoi(sdf[8].c_str())- snp;
        
           if(diff > snp)
           snp = diff;

           //int rs = (int)((double) (snp*100)/(double)(atoi(sdf[8].c_str())));  
 
           int rs = (int)((double) (snp*100)/(double)(atoi(sdf[8].c_str())+snp));

           if(m.count(rs) == 0)
             m[rs] = 1;
           else 
             m[rs]++; 
         }
     }
}

int main(int argc, char* argv[])
{
  ifstream in(argv[1]);
  char line[5000],recline[5000];
  while(in)
  {
     in.getline(line,5000);
     if(in)
     {
        readsnpfile1(line);
     }
  }
  in.close();

  for (map<int,int>:: iterator it = m.begin(); it != m.end(); it++ )
  cout << it->first<<"\t"<< it->second<< endl;
}
