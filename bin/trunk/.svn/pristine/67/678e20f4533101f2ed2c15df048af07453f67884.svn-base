#include "SXDBIndelgChecker.h" 


void  CSXDBIndelgChecker::ReverseComplementINDEL(char *seq, char *pcRevSeq)
{

     int j=0;
     for (int i = strlen(seq)-1; i >= 0; i--)
     {
         switch (seq[i])
         {
                 case 'A': pcRevSeq[j++] = 'T'; break;
                 case 'C': pcRevSeq[j++] = 'G'; break;
                 case 'G': pcRevSeq[j++] = 'C'; break;
                 case 'T': pcRevSeq[j++] = 'A';// break;
        }
   }
}


bool CSXDBIndelgChecker::matchINDEL(const char allele ,const char orient,const char* db)
{
     m_dbvec_INDEL.clear();

     if(db)
     {
        m_lendb_INDEL = strlen(db); m_ndb_INDEL = new char[m_lendb_INDEL+1]; strcpy(m_ndb_INDEL,db);        
        m_pch1_INDEL = strtok(m_ndb_INDEL,"/");

        while(m_pch1_INDEL != NULL)
        {
              if (orient=='-'){
                  char *pcRevSeq = new char[strlen(m_pch1_INDEL)+1];
                  ReverseComplementINDEL(m_pch1_INDEL,pcRevSeq); m_dbvec_INDEL.push_back(pcRevSeq);
                  delete[] pcRevSeq;          
              }
              else m_dbvec_INDEL.push_back(m_pch1_INDEL);
              
              m_pch1_INDEL = strtok(NULL,"/");
        }
        delete[] m_ndb_INDEL;
     }

     for (unsigned int ii =0; ii< m_dbvec_INDEL.size();ii++) {
          if (m_dbvec_INDEL[ii][0]==allele) return 1;      
     }
  
     return 0;
}


inline bool CSXDBIndelgChecker::MapAllele( Key* location ,const char variant,  FILE* pf)
{
     m_bFound_INDEL = false;
     m_InDelret = m_dbInDelMap.equal_range(location);

     for (m_dbInDelit=m_InDelret.first; m_dbInDelit!=m_InDelret.second; ++m_dbInDelit)
     {    
         if (matchINDEL(variant,(m_dbInDelit->second)->orient,(m_dbInDelit->second)->allele))
         {
             if (!m_bFound_INDEL){ m_bFound_INDEL = true; fprintf(pf,"%s",(m_dbInDelit->second)->rsid); }
             else { fprintf(pf,",%s",(m_dbInDelit->second)->rsid); }
         }
     }
     return m_bFound_INDEL;
}


bool CSXDBIndelgChecker::OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf)
{    
    m_bdbINDEL = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = m_dbInDelMap.find (snpkey);

    if (elefound != m_dbInDelMap.end ())
    {
        if (MapAllele(snpkey,cVariant,pf)) m_bdbINDEL = true;        
        else fprintf(pf,"-");         
    }
    else fprintf(pf,"-");
    
    return m_bdbINDEL;
}


inline bool CSXDBIndelgChecker::MapAllele( Key* location ,const char &variant, char acRsid[])
{
     m_bFound_INDEL = false; m_InDelret = m_dbInDelMap.equal_range(location);

     for (m_dbInDelit=m_InDelret.first; m_dbInDelit!=m_InDelret.second; ++m_dbInDelit)
     {
         if (matchINDEL(variant,(m_dbInDelit->second)->orient,(m_dbInDelit->second)->allele))
         {
             if (!m_bFound_INDEL){ m_bFound_INDEL = true; sprintf(acRsid,"%s",(m_dbInDelit->second)->rsid); }
             else { strcat(acRsid,","); strcat(acRsid,(m_dbInDelit->second)->rsid); }
         }
     }
     return m_bFound_INDEL;
}


bool CSXDBIndelgChecker::OutputINDELNovel(int chrno, size_t offset,const char cVariant, 
                                          FILE* pf, char acRsid[])
{
    m_bdbINDEL = false; Key * snpkey = new Key(chrno, offset);
    dbmap::const_iterator elefound = m_dbInDelMap.find (snpkey);

    if (elefound != m_dbInDelMap.end ()) {
        if (MapAllele(snpkey,cVariant,acRsid)) m_bdbINDEL = true;        
        else sprintf(acRsid,"-");        
    }
    else sprintf(acRsid,"-");    

    return m_bdbINDEL;
}


bool CSXDBIndelgChecker::GenerateDBINDELMap(FILE *dbfile)
{
    if (!dbfile) return false;
    char dbline[40960], *pch, *lpc[12]; int i; int chr_num;

    int nLen;

    while (fgets(dbline, sizeof(dbline), dbfile))
    {
        i=0; pch = strtok (dbline,"\t");

        while(pch != NULL){
          lpc[i++]=pch; if (i>12) break; pch = strtok (NULL,"\t");
        }

        if (lpc[0][0] == 'X') chr_num = 23;        
        else if (lpc[0][0] == 'Y') chr_num = 24;        
        else chr_num = atoi(lpc[0]);
       
        nLen = sizeof(lpc[4]);   

        if (nLen > 1)
        {
            for (int i = atoi(lpc[1]); i < atoi(lpc[2])+1; i++)
            {
                Key * key = new Key(chr_num, i);
                Value * value = new Value(lpc[8],lpc[3][0], lpc[5]);                  
                m_dbInDelMap.insert(pair<Key *, Value *>(key, value));
            }      
        }  
        else
        { 
            Key * key = new Key(chr_num, (size_t)atol(lpc[10])); // Key(int _chronum, size_t  _offset):chronum(_chronum),offset(_offset)
            Value * value = new Value(lpc[8],lpc[3][0], lpc[5]); // Value(const char* _rsid, const char * _orient, const char * _allele)
            m_dbInDelMap.insert(pair<Key *, Value *>(key, value));
        }

        memset(lpc,0,sizeof(lpc));
    }

    return true;
}


void CSXDBIndelgChecker::ClrDBINDELMap()
{
    m_dbInDelMap.clear();
}
