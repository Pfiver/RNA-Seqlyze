/* encodeRegionInfo.c was originally generated by the autoSql program, which also 
 * generated encodeRegionInfo.h and encodeRegionInfo.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "encode/encodeRegionInfo.h"


void encodeRegionInfoStaticLoad(char **row, struct encodeRegionInfo *ret)
/* Load a row from encodeRegionInfo table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->name = row[0];
ret->descr = row[1];
}

struct encodeRegionInfo *encodeRegionInfoLoad(char **row)
/* Load a encodeRegionInfo from row fetched with select * from encodeRegionInfo
 * from database.  Dispose of this with encodeRegionInfoFree(). */
{
struct encodeRegionInfo *ret;

AllocVar(ret);
ret->name = cloneString(row[0]);
ret->descr = cloneString(row[1]);
return ret;
}

struct encodeRegionInfo *encodeRegionInfoLoadAll(char *fileName) 
/* Load all encodeRegionInfo from a whitespace-separated file.
 * Dispose of this with encodeRegionInfoFreeList(). */
{
struct encodeRegionInfo *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileRow(lf, row))
    {
    el = encodeRegionInfoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct encodeRegionInfo *encodeRegionInfoLoadAllByChar(char *fileName, char chopper) 
/* Load all encodeRegionInfo from a chopper separated file.
 * Dispose of this with encodeRegionInfoFreeList(). */
{
struct encodeRegionInfo *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = encodeRegionInfoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct encodeRegionInfo *encodeRegionInfoCommaIn(char **pS, struct encodeRegionInfo *ret)
/* Create a encodeRegionInfo out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new encodeRegionInfo */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->descr = sqlStringComma(&s);
*pS = s;
return ret;
}

void encodeRegionInfoFree(struct encodeRegionInfo **pEl)
/* Free a single dynamically allocated encodeRegionInfo such as created
 * with encodeRegionInfoLoad(). */
{
struct encodeRegionInfo *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->descr);
freez(pEl);
}

void encodeRegionInfoFreeList(struct encodeRegionInfo **pList)
/* Free a list of dynamically allocated encodeRegionInfo's */
{
struct encodeRegionInfo *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    encodeRegionInfoFree(&el);
    }
*pList = NULL;
}

void encodeRegionInfoOutput(struct encodeRegionInfo *el, FILE *f, char sep, char lastSep) 
/* Print out encodeRegionInfo.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->descr);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

