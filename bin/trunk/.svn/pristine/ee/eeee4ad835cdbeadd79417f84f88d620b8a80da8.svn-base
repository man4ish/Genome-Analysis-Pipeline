#ifndef _SXAPPENDCYTPBAND_H
#define _SXAPPENDCYTPBAND_H

#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string.h>

    typedef struct _Cytoband
    {
         int nStart;
         int nStop;
         char *pcInfo;          

         _Cytoband()
         {
              nStart=nStop=0;
              pcInfo=NULL; 
         }

          
         ~_Cytoband()
         {
              if (pcInfo) free(pcInfo); 
         }

         inline void insert(int _nStart, int _nStop, char *_pcInfo)
         {
              this->nStart = _nStart;
              this->nStop = _nStop;
              this->pcInfo = (char*)malloc(strlen(_pcInfo)+1);
              strcpy(this->pcInfo,_pcInfo);
         } 
           
    } stCytoband;


    typedef struct _compareStartOffset{
        bool operator()(const stCytoband *ps1, const stCytoband *ps2)
        {
             if (ps1->nStart != ps2->nStart)
                 return ps1->nStart < ps2->nStart;
             else
                 return ps1->nStop < ps2->nStop;
        }
    };
 

    typedef std::vector<stCytoband*> CYTOBANDS_VEC;

    struct Cytoband_Chromosome
    {
         CYTOBANDS_VEC CytobandVec;
         CYTOBANDS_VEC::iterator CytobandItr;

         inline void sortCytobands()
         {
             sort(this->CytobandVec.begin(),this->CytobandVec.end(),_compareStartOffset());
         }

         inline void insertCytoband(int _nStart, int _nStop, char *_pcInfo)
         {
             stCytoband *pstCytoband = new stCytoband();
             pstCytoband->insert(_nStart,_nStop, _pcInfo);  
             this->CytobandVec.push_back(pstCytoband);          
         }               

         inline void printAll(const char *pcFile)   
         {
              FILE *pf = fopen(pcFile,"w");       

              for (this->CytobandItr = this->CytobandVec.begin(); this->CytobandItr != this->CytobandVec.end(); this->CytobandItr++) 
                   fprintf(pf,"%d\t%d\t%s\n",(*this->CytobandItr)->nStart,(*this->CytobandItr)->nStop,(*this->CytobandItr)->pcInfo);

              fclose(pf);
         }

         inline int search(int nOffset, char **lpcInfo, unsigned _maxCytobands)
         {       
              unsigned cnt = 0;

              CYTOBANDS_VEC::iterator start_itr = this->CytobandVec.begin();

              if (start_itr == this->CytobandVec.end()) return cnt;

              unsigned int unTotalRecs = this->CytobandVec.size();
              unsigned int unPrevStartPos = 0;
              unsigned int unStartPos = unTotalRecs/2;
              unsigned int unLastPos = unTotalRecs;
              unsigned int unMaxPlus1 = unTotalRecs+1;
              int nLoop=0;

              do {
                  start_itr = this->CytobandVec.begin()+unStartPos;
           
                  if (nOffset < (*start_itr)->nStart)
                  {
                      unLastPos = unStartPos;
                      unStartPos = (unPrevStartPos+unLastPos)/2;
                  }     
                  else if (nOffset > (*start_itr)->nStop) 
                  {
                     if (!(unLastPos >(unStartPos+1)))
                           break;

                     unPrevStartPos = unStartPos;
                     unStartPos=(unStartPos+unLastPos)/2;
                  }
                  else
                  {
                     break;
                  }
                  nLoop++;                 
              } while (unStartPos > 0 && unStartPos < unMaxPlus1 && nLoop < 100);

              if (!(nLoop < 100)){
                 fprintf(stderr,"offset = %u\n", nOffset);
                 start_itr = this->CytobandVec.begin();
              }

              if (unStartPos >(_maxCytobands-1)) start_itr -=_maxCytobands;//move _maxCytobands steps backwards to cater for offset repetition
              else start_itr = this->CytobandVec.begin();
    
              for (this->CytobandItr=start_itr; this->CytobandItr != this->CytobandVec.end(); this->CytobandItr++)  
              {
                   if (nOffset < (*this->CytobandItr)->nStart) break;

                   if (!(nOffset > (*this->CytobandItr)->nStop))
                   {
                       lpcInfo[cnt++] = (*this->CytobandItr)->pcInfo;
                   } 
              }
              return cnt;
         }
    }; //End of struct Cytoband_Chromosome

    
    typedef struct _compareChrNum{
         bool operator()(int s1, int s2)
         {
              return s1<s2;
         }
    };

    typedef std::map<int, Cytoband_Chromosome*, _compareChrNum> CYTOBANDS_MAP;
  
      
    FILE *g_pfIn;
    FILE *g_pfCytobands;
    FILE *g_pfOut; 

    CYTOBANDS_MAP g_Cytobands_Map;
    CYTOBANDS_MAP::iterator g_Cytobands_Itr;

#endif
