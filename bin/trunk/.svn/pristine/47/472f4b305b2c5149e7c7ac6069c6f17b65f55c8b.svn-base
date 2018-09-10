/* 
 * File:   SXDBIndelgChecker.h
 * Author: Hwah Shih Yiew
 *
 * Created on December 4, 2009, 11:01 AM
 */

#ifndef _SXDBINDELGCHECKER_H
#define	_SXDBINDELGCHECKER_H


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

struct Key
{
    int gi;
    size_t offset;
    Key(int  _gi, size_t  _offset):gi(_gi),offset(_offset) {}
};

struct Value
{
    char orient, *allele, *rsid;

    Value(const char* _rsid, const char _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1]; strcpy(rsid, _rsid);
        orient = _orient;
        allele = new char [strlen(_allele) + 1]; strcpy(allele, _allele);
    }

    ~Value()
    {
        if (rsid != NULL) delete [] rsid; if (allele != NULL) delete [] allele;
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi) return ( _lhs->gi < _rhs->gi);
        return  ( _lhs->offset < _rhs->offset );
    }
};


typedef multimap< Key*,Value*,ltstr> dbmap;

class CSXDBIndelgChecker
{
public:
    CSXDBIndelgChecker(){};
    virtual ~CSXDBIndelgChecker(){};

    bool GenerateDBINDELMap(FILE *dbfile);
    bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf);
    bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf,
                          char acRsid[]);
    void ClrDBINDELMap();


private:
    void ReverseComplementINDEL(char *seq, char *pcRevSeq);
    bool matchINDEL(const char allele ,const char orient,const char* db);
    inline bool MapAllele( Key* location ,const char variant,  FILE* pf);
    inline bool MapAllele( Key* location ,const char &variant, char acRsid[]);


private:
    dbmap m_dbInDelMap; 
    dbmap::iterator m_dbInDelit; 
    pair< dbmap ::iterator,dbmap::iterator> m_InDelret;
    bool m_bFound_INDEL, m_bdbINDEL;

    int m_lendb_INDEL; char *m_ndb_INDEL, *m_pch1_INDEL;
    vector <string>m_dbvec_INDEL;

public:
};
#endif	/* _SXDBINDELGCHECKER_H */

