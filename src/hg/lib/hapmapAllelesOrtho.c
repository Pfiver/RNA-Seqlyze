/* hapmapAllelesOrtho.c was originally generated by the autoSql program, which also 
 * generated hapmapAllelesOrtho.h and hapmapAllelesOrtho.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "hapmapAllelesOrtho.h"


void hapmapAllelesOrthoStaticLoad(char **row, struct hapmapAllelesOrtho *ret)
/* Load a row from hapmapAllelesOrtho table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->refUCSC = row[6];
ret->observed = row[7];
ret->orthoChrom = row[8];
ret->orthoStart = sqlUnsigned(row[9]);
ret->orthoEnd = sqlUnsigned(row[10]);
safecpy(ret->orthoStrand, sizeof(ret->orthoStrand), row[11]);
safecpy(ret->orthoAllele, sizeof(ret->orthoAllele), row[12]);
}

struct hapmapAllelesOrtho *hapmapAllelesOrthoLoad(char **row)
/* Load a hapmapAllelesOrtho from row fetched with select * from hapmapAllelesOrtho
 * from database.  Dispose of this with hapmapAllelesOrthoFree(). */
{
struct hapmapAllelesOrtho *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->refUCSC = cloneString(row[6]);
ret->observed = cloneString(row[7]);
ret->orthoChrom = cloneString(row[8]);
ret->orthoStart = sqlUnsigned(row[9]);
ret->orthoEnd = sqlUnsigned(row[10]);
safecpy(ret->orthoStrand, sizeof(ret->orthoStrand), row[11]);
safecpy(ret->orthoAllele, sizeof(ret->orthoAllele), row[12]);
return ret;
}

struct hapmapAllelesOrtho *hapmapAllelesOrthoLoadAll(char *fileName) 
/* Load all hapmapAllelesOrtho from a whitespace-separated file.
 * Dispose of this with hapmapAllelesOrthoFreeList(). */
{
struct hapmapAllelesOrtho *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[13];

while (lineFileRow(lf, row))
    {
    el = hapmapAllelesOrthoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct hapmapAllelesOrtho *hapmapAllelesOrthoLoadAllByChar(char *fileName, char chopper) 
/* Load all hapmapAllelesOrtho from a chopper separated file.
 * Dispose of this with hapmapAllelesOrthoFreeList(). */
{
struct hapmapAllelesOrtho *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[13];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = hapmapAllelesOrthoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct hapmapAllelesOrtho *hapmapAllelesOrthoCommaIn(char **pS, struct hapmapAllelesOrtho *ret)
/* Create a hapmapAllelesOrtho out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new hapmapAllelesOrtho */
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
ret->refUCSC = sqlStringComma(&s);
ret->observed = sqlStringComma(&s);
ret->orthoChrom = sqlStringComma(&s);
ret->orthoStart = sqlUnsignedComma(&s);
ret->orthoEnd = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->orthoStrand, sizeof(ret->orthoStrand));
sqlFixedStringComma(&s, ret->orthoAllele, sizeof(ret->orthoAllele));
*pS = s;
return ret;
}

void hapmapAllelesOrthoFree(struct hapmapAllelesOrtho **pEl)
/* Free a single dynamically allocated hapmapAllelesOrtho such as created
 * with hapmapAllelesOrthoLoad(). */
{
struct hapmapAllelesOrtho *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->refUCSC);
freeMem(el->observed);
freeMem(el->orthoChrom);
freez(pEl);
}

void hapmapAllelesOrthoFreeList(struct hapmapAllelesOrtho **pList)
/* Free a list of dynamically allocated hapmapAllelesOrtho's */
{
struct hapmapAllelesOrtho *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    hapmapAllelesOrthoFree(&el);
    }
*pList = NULL;
}

void hapmapAllelesOrthoOutput(struct hapmapAllelesOrtho *el, FILE *f, char sep, char lastSep) 
/* Print out hapmapAllelesOrtho.  Separate fields with sep. Follow last field with lastSep. */
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
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->refUCSC);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->observed);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->orthoChrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->orthoStart);
fputc(sep,f);
fprintf(f, "%u", el->orthoEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->orthoStrand);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->orthoAllele);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

