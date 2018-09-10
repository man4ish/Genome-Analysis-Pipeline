/* 
 * File:   SXDBSNPgChecker.h
 * Author: Hwah Shih Yiew
 *
 * Created on December 3, 2009, 2:06 PM
 */

#ifndef _SXDBSNPGCHECKER_H
#define	_SXDBSNPGCHECKER_H

bool GenerateDBSNPMap(FILE *dbfile);
bool GenerateDBSNPMap_OdbSnpfmt(FILE *dbfile);
bool OutputSNPNovel(int chrno, size_t offset,const char cVariant, FILE* pf);
bool OutputSNPNovel2(int chrno, size_t offset,const char cVariant, FILE* pf);
bool OutputSNPNovel(int chrno, size_t offset,const char cVariant, FILE* pf,
                    char acRsid[]);
void GetSNPStatistics(const char* flname);
void ClrDBSNPMap();
#endif	/* _SXDBSNPGCHECKER_H */

