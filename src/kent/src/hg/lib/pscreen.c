/* pscreen.c was originally generated by the autoSql program, which also 
 * generated pscreen.h and pscreen.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "pscreen.h"


struct pscreen *pscreenLoad(char **row)
/* Load a pscreen from row fetched with select * from pscreen
 * from database.  Dispose of this with pscreenFree(). */
{
struct pscreen *ret;

AllocVar(ret);
ret->geneCount = sqlUnsigned(row[7]);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
strcpy(ret->strand, row[5]);
ret->stockNumber = sqlUnsigned(row[6]);
{
int sizeOne;
sqlStringDynamicArray(row[8], &ret->geneIds, &sizeOne);
assert(sizeOne == ret->geneCount);
}
return ret;
}

struct pscreen *pscreenLoadAll(char *fileName) 
/* Load all pscreen from a whitespace-separated file.
 * Dispose of this with pscreenFreeList(). */
{
struct pscreen *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[9];

while (lineFileRow(lf, row))
    {
    el = pscreenLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct pscreen *pscreenLoadAllByChar(char *fileName, char chopper) 
/* Load all pscreen from a chopper separated file.
 * Dispose of this with pscreenFreeList(). */
{
struct pscreen *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[9];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = pscreenLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct pscreen *pscreenCommaIn(char **pS, struct pscreen *ret)
/* Create a pscreen out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new pscreen */
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
ret->stockNumber = sqlUnsignedComma(&s);
ret->geneCount = sqlUnsignedComma(&s);
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->geneIds, ret->geneCount);
for (i=0; i<ret->geneCount; ++i)
    {
    ret->geneIds[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
*pS = s;
return ret;
}

void pscreenFree(struct pscreen **pEl)
/* Free a single dynamically allocated pscreen such as created
 * with pscreenLoad(). */
{
struct pscreen *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
/* All strings in geneIds are allocated at once, so only need to free first. */
if (el->geneIds != NULL)
    freeMem(el->geneIds[0]);
freeMem(el->geneIds);
freez(pEl);
}

void pscreenFreeList(struct pscreen **pList)
/* Free a list of dynamically allocated pscreen's */
{
struct pscreen *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    pscreenFree(&el);
    }
*pList = NULL;
}

void pscreenOutput(struct pscreen *el, FILE *f, char sep, char lastSep) 
/* Print out pscreen.  Separate fields with sep. Follow last field with lastSep. */
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
fputc(sep,f);
fprintf(f, "%u", el->stockNumber);
fputc(sep,f);
fprintf(f, "%u", el->geneCount);
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->geneCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->geneIds[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

