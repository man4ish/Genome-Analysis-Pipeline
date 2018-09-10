/* 
 * File:   SXDBIndelgChecker.h
 * Author: Hwah Shih Yiew
 *
 * Created on December 4, 2009, 11:01 AM
 */

#ifndef _SXDBINDELGCHECKER_H
#define	_SXDBINDELGCHECKER_H

bool GenerateDBINDELMap(FILE *dbfile);
bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf);
bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf);
bool OutputINDELNovel(int chrno, size_t offset,const char cVariant, FILE* pf,
                      char acRsid[]);
void ClrDBINDELMap();

#endif	/* _SXDBINDELGCHECKER_H */

