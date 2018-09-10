#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
map<int,int> m;
float cutoff;
void manual()
{
     cout<< "\n\nSXFiltermean3sd\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXFiltermean3sd <cutoff> <Input File>  \n\n";
}

void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chr"))
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
           int chrnum;
           if(sdf[0]==  "X")
           chrnum = 23;
           else if (sdf[0] ==  "Y")
           chrnum = 24;
           else 
           chrnum = atoi(sdf[0].c_str());
           if( atoi(sdf[8].c_str()) <= cutoff )
            cout << line1<< endl;
         }
     }
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    manual();
    exit(1);
  }
  
  cutoff = atof(argv[1]);

  ifstream in(argv[2]);
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
}
