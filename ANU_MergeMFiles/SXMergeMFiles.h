#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <map>

typedef struct _stKey
{
    int nOffset;
    char *pcAllele;

    _stKey(){
        nOffset = 0; pcAllele=NULL;  
    }

    ~_stKey(){
        if (pcAllele) delete[] pcAllele;
    }
}stKey;


typedef struct _comparekey
{
    bool operator()(const stKey *p1, const stKey *p2)const
    {
         if (p1->nOffset != p2->nOffset)
            return p1->nOffset < p2->nOffset;

         return strcmp(p1->pcAllele, p2->pcAllele) < 0;  
    }
};


typedef std::map<const stKey*, int, _comparekey> M_List;
