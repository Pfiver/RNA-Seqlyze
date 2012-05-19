/* knownCanonical.c was originally generated by the autoSql program, which also 
 * generated knownCanonical.h and knownCanonical.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "knownCanonical.h"


void knownCanonicalStaticLoad(char **row, struct knownCanonical *ret)
/* Load a row from knownCanonical table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlSigned(row[1]);
ret->chromEnd = sqlSigned(row[2]);
ret->clusterId = sqlSigned(row[3]);
ret->transcript = row[4];
ret->protein = row[5];
}

struct knownCanonical *knownCanonicalLoad(char **row)
/* Load a knownCanonical from row fetched with select * from knownCanonical
 * from database.  Dispose of this with knownCanonicalFree(). */
{
struct knownCanonical *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlSigned(row[1]);
ret->chromEnd = sqlSigned(row[2]);
ret->clusterId = sqlSigned(row[3]);
ret->transcript = cloneString(row[4]);
ret->protein = cloneString(row[5]);
return ret;
}

struct knownCanonical *knownCanonicalLoadAll(char *fileName) 
/* Load all knownCanonical from a whitespace-separated file.
 * Dispose of this with knownCanonicalFreeList(). */
{
struct knownCanonical *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[6];

while (lineFileRow(lf, row))
    {
    el = knownCanonicalLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct knownCanonical *knownCanonicalLoadAllByChar(char *fileName, char chopper) 
/* Load all knownCanonical from a chopper separated file.
 * Dispose of this with knownCanonicalFreeList(). */
{
struct knownCanonical *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[6];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = knownCanonicalLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct knownCanonical *knownCanonicalCommaIn(char **pS, struct knownCanonical *ret)
/* Create a knownCanonical out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new knownCanonical */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlSignedComma(&s);
ret->chromEnd = sqlSignedComma(&s);
ret->clusterId = sqlSignedComma(&s);
ret->transcript = sqlStringComma(&s);
ret->protein = sqlStringComma(&s);
*pS = s;
return ret;
}

void knownCanonicalFree(struct knownCanonical **pEl)
/* Free a single dynamically allocated knownCanonical such as created
 * with knownCanonicalLoad(). */
{
struct knownCanonical *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->transcript);
freeMem(el->protein);
freez(pEl);
}

void knownCanonicalFreeList(struct knownCanonical **pList)
/* Free a list of dynamically allocated knownCanonical's */
{
struct knownCanonical *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    knownCanonicalFree(&el);
    }
*pList = NULL;
}

void knownCanonicalOutput(struct knownCanonical *el, FILE *f, char sep, char lastSep) 
/* Print out knownCanonical.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->chromStart);
fputc(sep,f);
fprintf(f, "%d", el->chromEnd);
fputc(sep,f);
fprintf(f, "%d", el->clusterId);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->transcript);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->protein);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

