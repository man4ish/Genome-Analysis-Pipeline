#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <vector>

using namespace std;
FILE *g_pfMaskDens[24];

unsigned int GetRepeatregion(int i, unsigned int offset)
{
    offset--;
    unsigned int g_unReads; unsigned char g_acBuf[2];   

    fseek(g_pfMaskDens[i],2*offset,SEEK_SET);   
    fread(g_acBuf,2,1,g_pfMaskDens[i]);

    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|=g_acBuf[0];
    
    return g_unReads;
}
 
int main(int argc, char* argv[])
{
   if (argc < 4)
   {
        fprintf(stdout,"Parameters not sufficient...\n");
        fprintf(stdout,"Usage: %s\n",argv[0]);
        fprintf(stdout,"\t<Input filename> <Masked Density files path> <Output filename>\n");
        exit(9);
   }

   char acFileName[1024];

   ifstream in(argv[1]);
   
   for (int i=0; i < 24; i++){
       sprintf(acFileName,"%s/Chr%02d_mask_den",argv[2],i+1);
       g_pfMaskDens[i]=fopen(acFileName,"rb");
       if (!g_pfMaskDens[i]) {fprintf(stdout,"Failed to open %s ...\n",acFileName); exit(9);} 
   }

   FILE *pfOut=fopen(argv[3],"w"); 
   if (!pfOut) {fprintf(stdout,"Failed to open %s ...\n",argv[3]); exit(9);}

   char line[5000], recline[5000], *pch;
   int nTest1=0, nTest2=0, rpc, nIdx;

   while(in)
   {
      if (++nTest1 > 100000)
      {
         fprintf(stderr,"Total Recs Processed = %u x 0.1M\n",++nTest2);
         nTest1=0;
      }

      in.getline(line,5000); strcpy(recline,line);

      if(in)
      {        
         pch = strtok(line,"\t"); vector<const char*> vec;
         while(pch != NULL)
         {   
             vec.push_back(pch); pch = strtok(NULL,"\t"); 
         }
         
         if (vec[0][0]=='X') nIdx=22; 
         else if (vec[0][0]=='Y') nIdx=23;
         else nIdx = atoi(vec[0]) -1;
           
         if (GetRepeatregion(nIdx,atoi(vec[1]))==0) fprintf(pfOut,"%s\n",recline);         
     }  
      
   }


   for (int i=0; i < 24; i++){
       if (g_pfMaskDens[i]) fclose(g_pfMaskDens[i]);
   }


  return 0;
}
