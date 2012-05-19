/* fosEndPairs.c was originally generated by the autoSql program, which also 
 * generated fosEndPairs.h and fosEndPairs.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "fosEndPairs.h"


struct fosEndPairs *fosEndPairsLoad(char **row)
/* Load a fosEndPairs from row fetched with select * from fosEndPairs
 * from database.  Dispose of this with fosEndPairsFree(). */
{
struct fosEndPairs *ret;
int sizeOne,i;
char *s;

AllocVar(ret);
ret->lfCount = sqlUnsigned(row[8]);
ret->bin = sqlSigned(row[0]);
ret->chrom = cloneString(row[1]);
ret->chromStart = sqlUnsigned(row[2]);
ret->chromEnd = sqlUnsigned(row[3]);
ret->name = cloneString(row[4]);
ret->score = sqlUnsigned(row[5]);
strcpy(ret->strand, row[6]);
ret->pslTable = cloneString(row[7]);
sqlUnsignedDynamicArray(row[9], &ret->lfStarts, &sizeOne);
assert(sizeOne == ret->lfCount);
sqlUnsignedDynamicArray(row[10], &ret->lfSizes, &sizeOne);
assert(sizeOne == ret->lfCount);
sqlStringDynamicArray(row[11], &ret->lfNames, &sizeOne);
assert(sizeOne == ret->lfCount);
return ret;
}

struct fosEndPairs *fosEndPairsLoadAll(char *fileName) 
/* Load all fosEndPairs from a tab-separated file.
 * Dispose of this with fosEndPairsFreeList(). */
{
struct fosEndPairs *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[12];

while (lineFileRow(lf, row))
    {
    el = fosEndPairsLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct fosEndPairs *fosEndPairsLoadWhere(struct sqlConnection *conn, char *table, char *where)
/* Load all fosEndPairs from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with fosEndPairsFreeList(). */
{
struct fosEndPairs *list = NULL, *el;
struct dyString *query = dyStringNew(256);
struct sqlResult *sr;
char **row;

dyStringPrintf(query, "select * from %s", table);
if (where != NULL)
    dyStringPrintf(query, " where %s", where);
sr = sqlGetResult(conn, query->string);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = fosEndPairsLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
dyStringFree(&query);
return list;
}

struct fosEndPairs *fosEndPairsCommaIn(char **pS, struct fosEndPairs *ret)
/* Create a fosEndPairs out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new fosEndPairs */
{
char *s = *pS;
int i;

if (ret == NULL)
    AllocVar(ret);
ret->bin = sqlSignedComma(&s);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
ret->pslTable = sqlStringComma(&s);
ret->lfCount = sqlUnsignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->lfStarts, ret->lfCount);
for (i=0; i<ret->lfCount; ++i)
    {
    ret->lfStarts[i] = sqlUnsignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
s = sqlEatChar(s, '{');
AllocArray(ret->lfSizes, ret->lfCount);
for (i=0; i<ret->lfCount; ++i)
    {
    ret->lfSizes[i] = sqlUnsignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
s = sqlEatChar(s, '{');
AllocArray(ret->lfNames, ret->lfCount);
for (i=0; i<ret->lfCount; ++i)
    {
    ret->lfNames[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
*pS = s;
return ret;
}

void fosEndPairsFree(struct fosEndPairs **pEl)
/* Free a single dynamically allocated fosEndPairs such as created
 * with fosEndPairsLoad(). */
{
struct fosEndPairs *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->pslTable);
freeMem(el->lfStarts);
freeMem(el->lfSizes);
/* All strings in lfNames are allocated at once, so only need to free first. */
if (el->lfNames != NULL)
    freeMem(el->lfNames[0]);
freeMem(el->lfNames);
freez(pEl);
}

void fosEndPairsFreeList(struct fosEndPairs **pList)
/* Free a list of dynamically allocated fosEndPairs's */
{
struct fosEndPairs *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    fosEndPairsFree(&el);
    }
*pList = NULL;
}

void fosEndPairsOutput(struct fosEndPairs *el, FILE *f, char sep, char lastSep) 
/* Print out fosEndPairs.  Separate fields with sep. Follow last field with lastSep. */
{
int i;
fprintf(f, "%d", el->bin);
fputc(sep,f);
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
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->pslTable);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->lfCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->lfCount; ++i)
    {
    fprintf(f, "%u", el->lfStarts[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->lfCount; ++i)
    {
    fprintf(f, "%u", el->lfSizes[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->lfCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->lfNames[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(lastSep,f);
}

