#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
#include <map>

typedef struct _stKey
{
     int nChro;
     int nOffset;
     char acAllele[4];

     _stKey()
     {
          nChro=nOffset=0;
          acAllele[0]=0;           
     }
}stKey;

typedef struct _stData
{
    char *pcLeftSeq;
    char *pcRightSeq;

    _stData()
    {
         pcLeftSeq=pcRightSeq=NULL;
    } 

    ~_stData()
    {
        if (pcLeftSeq) free(pcLeftSeq);
        if (pcRightSeq) free(pcRightSeq);
    }
}stData;

struct _compareKey
{
    bool operator()(const stKey *p1, const stKey *p2)
    {
         if (p1->nChro != p2->nChro) return p1->nChro < p2->nChro;
         if (p1->nOffset != p2->nOffset) return p1->nOffset < p2->nOffset; 
         if (p1->acAllele[0] != p2->acAllele[0]) return p1->acAllele[0] < p2->acAllele[0];
         return p1->acAllele[2] < p2->acAllele[2];
    }
};

typedef std::map<stKey*, stData*, _compareKey> VSEQ_MAP;

VSEQ_MAP g_VSEQ_Map;
VSEQ_MAP::iterator g_itr;
*/
FILE *g_pfIn, *g_pfVSeq, *g_pfOut;


static void banner(char *argv[])
{
    printf("Synamatix SXChkAppendVSeq Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: \n", argv[0]);
    printf("\t<Input filename> <Var. Sequence filename> <Output filename>\n\n");
}


void ProcessSeqRec()
{
    char acbuf[40960], acVSeq[4096], *pChr=NULL,*pcLast=NULL; 
    //stKey *pstKey; stData *pstData;
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfIn)){

       if(!fgets(acbuf,sizeof(acbuf),g_pfIn)) break;  

       pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
       pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

       if (acbuf[0]==0)  fprintf(g_pfOut,"\n");
       else if (acbuf[0]=='#') fprintf(g_pfOut,"%s\n",acbuf);
       else if (acbuf[0]=='C') { fprintf(g_pfOut,"%s\tLeft\tRight\n",acbuf); break;}
       else { fseek(g_pfIn,0,SEEK_SET); break;}
    }

    while(!feof(g_pfVSeq))
    { 
         if(!fgets(acbuf,sizeof(acbuf),g_pfVSeq)) break;

          pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
          pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

          if (acbuf[0]==0||acbuf[0]=='#') continue;
          else if (acbuf[0]=='C') break;   
          else { fseek(g_pfVSeq,0,SEEK_SET); break;}         
    }

    while(!feof(g_pfIn) && !feof(g_pfVSeq))
    {
          if (++ulTest1 == 100000)
          {
              fprintf(stdout,"Total Recs Processed = %u x 0.1M\n", ++ulTest2);
              ulTest1=0;
          }

          if(!fgets(acbuf,sizeof(acbuf),g_pfIn)) break; 

          pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
          pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

          fprintf(g_pfOut,"%s\t",acbuf); 
            
          if(!fgets(acbuf,sizeof(acbuf),g_pfVSeq)) break; 

          pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
          pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

          pChr = strtok(acbuf,"\t"); //Chro
          pChr = strtok(NULL,"\t"); //Start  
          pChr = strtok(NULL,"\t"); //End

          pChr = strtok(NULL,"\t"); //Left
          fprintf(g_pfOut,"%s\t",pChr);

          pChr = strtok(NULL,"\t"); //Right
          fprintf(g_pfOut,"%s\n",pChr);

          /*
          pstKey = new stKey; pstData=new stData;

          pChr = strtok(acbuf,"\t"); pstKey->nChro = atoi(pChr);
          pChr = strtok(NULL,"\t");  pstKey->nOffset = atoi(pChr);
          pChr = strtok(NULL,"\t");
          pChr = strtok(NULL,"\t");   
          pChr = strtok(NULL,"\t");
          */ 
    }
}



int main(int argc, char *argv[])
{
    if (argc != 4){ printf("Invalid parameters\n");banner(argv); exit(9);}

    g_pfIn = fopen(argv[1],"r"); if (!g_pfIn) {printf("Failed to open %s ...\n",argv[1]); goto ExitRtn;}
    g_pfVSeq = fopen(argv[2],"r"); if (!g_pfVSeq){printf("Failed to open %s ...\n",argv[2]); goto ExitRtn;}
    g_pfOut = fopen(argv[3],"w"); if (!g_pfOut) {printf("Failed to open %s ...\n",argv[3]); goto ExitRtn;}

    printf("Start Appending VSeq Data ...\n"); ProcessSeqRec();
    printf("Done ...\n");

ExitRtn:
    ;


}





