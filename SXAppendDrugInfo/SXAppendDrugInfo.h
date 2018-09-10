#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <map>

typedef struct _stDrugInfo{
    char *pcEvidence;
    char *pcAnnotation;
    char *pcDrugs;
    char *pcDrugClasses;
    char *pcDisease; 

    _stDrugInfo(){
       pcEvidence = NULL;
       pcAnnotation = NULL;
       pcDrugs = NULL;
       pcDrugClasses = NULL;
       pcDisease = NULL;    
    }

    ~_stDrugInfo(){
       if (pcEvidence) delete[] pcEvidence;
       if (pcAnnotation) delete[] pcAnnotation;
       if (pcDrugs) delete[] pcDrugs;
       if (pcDrugClasses) delete[] pcDrugClasses;
       if (pcDisease) delete [] pcDisease;
    }
} stDrugInfo;


struct _comparekey
{
    bool operator()(const char *p1, const char *p2) const
    {
         return strcmp(p1,p2) <  0;
    }
};


typedef std::map<const char*, stDrugInfo*, _comparekey> DrugInfo_List;
