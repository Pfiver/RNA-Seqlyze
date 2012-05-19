/* syntenySanger.c was originally generated by the autoSql program, which also 
 * generated syntenySanger.h and syntenySanger.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "syntenySanger.h"


void syntenySangerStaticLoad(char **row, struct syntenySanger *ret)
/* Load a row from syntenySanger table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
strcpy(ret->strand, row[5]);
}

struct syntenySanger *syntenySangerLoad(char **row)
/* Load a syntenySanger from row fetched with select * from syntenySanger
 * from database.  Dispose of this with syntenySangerFree(). */
{
struct syntenySanger *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
strcpy(ret->strand, row[5]);
return ret;
}

struct syntenySanger *syntenySangerLoadAll(char *fileName) 
/* Load all syntenySanger from a tab-separated file.
 * Dispose of this with syntenySangerFreeList(). */
{
struct syntenySanger *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[6];

while (lineFileRow(lf, row))
    {
    el = syntenySangerLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct syntenySanger *syntenySangerLoadWhere(struct sqlConnection *conn, char *table, char *where)
/* Load all syntenySanger from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with syntenySangerFreeList(). */
{
struct syntenySanger *list = NULL, *el;
struct dyString *query = dyStringNew(256);
struct sqlResult *sr;
char **row;

dyStringPrintf(query, "select * from %s", table);
if (where != NULL)
    dyStringPrintf(query, " where %s", where);
sr = sqlGetResult(conn, query->string);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = syntenySangerLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
dyStringFree(&query);
return list;
}

struct syntenySanger *syntenySangerCommaIn(char **pS, struct syntenySanger *ret)
/* Create a syntenySanger out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new syntenySanger */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
*pS = s;
return ret;
}

void syntenySangerFree(struct syntenySanger **pEl)
/* Free a single dynamically allocated syntenySanger such as created
 * with syntenySangerLoad(). */
{
struct syntenySanger *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void syntenySangerFreeList(struct syntenySanger **pList)
/* Free a list of dynamically allocated syntenySanger's */
{
struct syntenySanger *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    syntenySangerFree(&el);
    }
*pList = NULL;
}

void syntenySangerOutput(struct syntenySanger *el, FILE *f, char sep, char lastSep) 
/* Print out syntenySanger.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->score);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

