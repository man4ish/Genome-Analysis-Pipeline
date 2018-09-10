#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <vector>
using namespace std;
unsigned int GetRepeatregion(int i, unsigned int offset)
{
    string filename;
    filename = "/archive/project/oneclick/syhwah/SXFilterRepeat/Reptmask_20100515/Chr"; //"Reptmask/Chr"; 
    std::stringstream out;
    out << i;       
    string s = filename;
    s  += out.str();
    s += "_mask_den";
    FILE* pf = fopen(s.c_str(),"rb");
    if (!pf) return 0;
    unsigned int g_unReads;   
    unsigned char g_acBuf[2];   
    fseek(pf,2*offset,SEEK_SET);   
    fread(g_acBuf,2,1,pf);
    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|=g_acBuf[0];
    fclose(pf);
    return g_unReads;
}
 
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        fprintf(stdout,"Parameters not sufficient...\n");
        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename>\n");
        exit(9);
    }

   ifstream in(argv[1]);
   char line[5000], recline[5000];
   while(in)
   {
      in.getline(line,5000);
      strcpy(recline,line);
      if(in)
      {
         char* pch;
         pch = strtok(line,"\t");
         vector<const char*> vec;
         while(pch != NULL)
         {   
             vec.push_back(pch);
             pch = strtok(NULL,"\t"); 
         }
         int rpc = GetRepeatregion(atoi(vec[0]), atoi(vec[1]));
         if(rpc == 0)
           cout  << recline << endl;
         
      }
   }
  return 0;
}
