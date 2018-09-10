#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

typedef struct _stDenFile{
     int fd;
     size_t FileSize;        
     char *pcData; 

     _stDenFile()
     {
         fd=0; FileSize=0; pcData=NULL;    
     }
} stDenFile;



FILE *g_pfIn=NULL, *g_pfOut=NULL;//,**g_apfMaskedDen=NULL,**g_apfReadDen=NULL,**g_apfReptDen=NULL;

stDenFile g_stReadDen[24];
stDenFile g_stReptDen[24];
//stDenFile g_stMaskedDen[24];

char *g_pcSampleID, g_cSNVType; float g_fCutOff,g_fMeans;
int g_nHetMinSR, g_nHomMinSR;

