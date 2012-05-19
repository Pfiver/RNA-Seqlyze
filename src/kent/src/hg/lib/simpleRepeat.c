/* simpleRepeat.c was originally generated by the autoSql program, which also 
 * generated simpleRepeat.h and simpleRepeat.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "simpleRepeat.h"


void simpleRepeatStaticLoad(char **row, struct simpleRepeat *ret)
/* Load a row from simpleRepeat table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->period = sqlUnsigned(row[4]);
ret->copyNum = atof(row[5]);
ret->consensusSize = sqlUnsigned(row[6]);
ret->perMatch = sqlUnsigned(row[7]);
ret->perIndel = sqlUnsigned(row[8]);
ret->score = sqlUnsigned(row[9]);
ret->A = sqlUnsigned(row[10]);
ret->C = sqlUnsigned(row[11]);
ret->G = sqlUnsigned(row[12]);
ret->T = sqlUnsigned(row[13]);
ret->entropy = atof(row[14]);
ret->sequence = row[15];
}

struct simpleRepeat *simpleRepeatLoad(char **row)
/* Load a simpleRepeat from row fetched with select * from simpleRepeat
 * from database.  Dispose of this with simpleRepeatFree(). */
{
struct simpleRepeat *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->period = sqlUnsigned(row[4]);
ret->copyNum = atof(row[5]);
ret->consensusSize = sqlUnsigned(row[6]);
ret->perMatch = sqlUnsigned(row[7]);
ret->perIndel = sqlUnsigned(row[8]);
ret->score = sqlUnsigned(row[9]);
ret->A = sqlUnsigned(row[10]);
ret->C = sqlUnsigned(row[11]);
ret->G = sqlUnsigned(row[12]);
ret->T = sqlUnsigned(row[13]);
ret->entropy = atof(row[14]);
ret->sequence = cloneString(row[15]);
return ret;
}

struct simpleRepeat *simpleRepeatLoadAll(char *fileName) 
/* Load all simpleRepeat from a tab-separated file.
 * Dispose of this with simpleRepeatFreeList(). */
{
struct simpleRepeat *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[16];

while (lineFileRow(lf, row))
    {
    el = simpleRepeatLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct simpleRepeat *simpleRepeatLoadWhere(struct sqlConnection *conn, char *table, char *where)
/* Load all simpleRepeat from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with simpleRepeatFreeList(). */
{
struct simpleRepeat *list = NULL, *el;
struct dyString *query = dyStringNew(256);
struct sqlResult *sr;
char **row;

dyStringPrintf(query, "select * from %s", table);
if (where != NULL)
    dyStringPrintf(query, " where %s", where);
sr = sqlGetResult(conn, query->string);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = simpleRepeatLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
dyStringFree(&query);
return list;
}

struct simpleRepeat *simpleRepeatCommaIn(char **pS, struct simpleRepeat *ret)
/* Create a simpleRepeat out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new simpleRepeat */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->period = sqlUnsignedComma(&s);
ret->copyNum = sqlFloatComma(&s);
ret->consensusSize = sqlUnsignedComma(&s);
ret->perMatch = sqlUnsignedComma(&s);
ret->perIndel = sqlUnsignedComma(&s);
ret->score = sqlUnsignedComma(&s);
ret->A = sqlUnsignedComma(&s);
ret->C = sqlUnsignedComma(&s);
ret->G = sqlUnsignedComma(&s);
ret->T = sqlUnsignedComma(&s);
ret->entropy = sqlFloatComma(&s);
ret->sequence = sqlStringComma(&s);
*pS = s;
return ret;
}

void simpleRepeatFree(struct simpleRepeat **pEl)
/* Free a single dynamically allocated simpleRepeat such as created
 * with simpleRepeatLoad(). */
{
struct simpleRepeat *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->sequence);
freez(pEl);
}

void simpleRepeatFreeList(struct simpleRepeat **pList)
/* Free a list of dynamically allocated simpleRepeat's */
{
struct simpleRepeat *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    simpleRepeatFree(&el);
    }
*pList = NULL;
}

void simpleRepeatOutput(struct simpleRepeat *el, FILE *f, char sep, char lastSep) 
/* Print out simpleRepeat.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%u", el->period);
fputc(sep,f);
fprintf(f, "%f", el->copyNum);
fputc(sep,f);
fprintf(f, "%u", el->consensusSize);
fputc(sep,f);
fprintf(f, "%u", el->perMatch);
fputc(sep,f);
fprintf(f, "%u", el->perIndel);
fputc(sep,f);
fprintf(f, "%u", el->score);
fputc(sep,f);
fprintf(f, "%u", el->A);
fputc(sep,f);
fprintf(f, "%u", el->C);
fputc(sep,f);
fprintf(f, "%u", el->G);
fputc(sep,f);
fprintf(f, "%u", el->T);
fputc(sep,f);
fprintf(f, "%f", el->entropy);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->sequence);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

