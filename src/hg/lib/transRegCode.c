/* transRegCode.c was originally generated by the autoSql program, which also 
 * generated transRegCode.h and transRegCode.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "transRegCode.h"


void transRegCodeStaticLoad(char **row, struct transRegCode *ret)
/* Load a row from transRegCode table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
ret->chipEvidence = row[5];
ret->consSpecies = sqlUnsigned(row[6]);
}

struct transRegCode *transRegCodeLoad(char **row)
/* Load a transRegCode from row fetched with select * from transRegCode
 * from database.  Dispose of this with transRegCodeFree(). */
{
struct transRegCode *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
ret->chipEvidence = cloneString(row[5]);
ret->consSpecies = sqlUnsigned(row[6]);
return ret;
}

struct transRegCode *transRegCodeLoadAll(char *fileName) 
/* Load all transRegCode from a whitespace-separated file.
 * Dispose of this with transRegCodeFreeList(). */
{
struct transRegCode *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[7];

while (lineFileRow(lf, row))
    {
    el = transRegCodeLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct transRegCode *transRegCodeLoadAllByChar(char *fileName, char chopper) 
/* Load all transRegCode from a chopper separated file.
 * Dispose of this with transRegCodeFreeList(). */
{
struct transRegCode *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[7];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = transRegCodeLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct transRegCode *transRegCodeCommaIn(char **pS, struct transRegCode *ret)
/* Create a transRegCode out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new transRegCode */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
ret->chipEvidence = sqlStringComma(&s);
ret->consSpecies = sqlUnsignedComma(&s);
*pS = s;
return ret;
}

void transRegCodeFree(struct transRegCode **pEl)
/* Free a single dynamically allocated transRegCode such as created
 * with transRegCodeLoad(). */
{
struct transRegCode *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->chipEvidence);
freez(pEl);
}

void transRegCodeFreeList(struct transRegCode **pList)
/* Free a list of dynamically allocated transRegCode's */
{
struct transRegCode *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    transRegCodeFree(&el);
    }
*pList = NULL;
}

void transRegCodeOutput(struct transRegCode *el, FILE *f, char sep, char lastSep) 
/* Print out transRegCode.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%s", el->chipEvidence);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->consSpecies);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

