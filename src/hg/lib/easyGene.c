/* easyGene.c was originally generated by the autoSql program, which also 
 * generated easyGene.h and easyGene.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "easyGene.h"


void easyGeneStaticLoad(char **row, struct easyGene *ret)
/* Load a row from easyGene table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
strcpy(ret->strand, row[5]);
ret->feat = row[6];
ret->R = atof(row[7]);
ret->frame = sqlUnsigned(row[8]);
ret->orf = row[9];
ret->startCodon = row[10];
ret->logOdds = atof(row[11]);
ret->descriptor = row[12];
ret->swissProt = row[13];
strcpy(ret->genbank, row[14]);
}

struct easyGene *easyGeneLoad(char **row)
/* Load a easyGene from row fetched with select * from easyGene
 * from database.  Dispose of this with easyGeneFree(). */
{
struct easyGene *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
strcpy(ret->strand, row[5]);
ret->feat = cloneString(row[6]);
ret->R = atof(row[7]);
ret->frame = sqlUnsigned(row[8]);
ret->orf = cloneString(row[9]);
ret->startCodon = cloneString(row[10]);
ret->logOdds = atof(row[11]);
ret->descriptor = cloneString(row[12]);
ret->swissProt = cloneString(row[13]);
strcpy(ret->genbank, row[14]);
return ret;
}

struct easyGene *easyGeneLoadAll(char *fileName) 
/* Load all easyGene from a whitespace-separated file.
 * Dispose of this with easyGeneFreeList(). */
{
struct easyGene *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileRow(lf, row))
    {
    el = easyGeneLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct easyGene *easyGeneLoadAllByChar(char *fileName, char chopper) 
/* Load all easyGene from a chopper separated file.
 * Dispose of this with easyGeneFreeList(). */
{
struct easyGene *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[15];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = easyGeneLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct easyGene *easyGeneCommaIn(char **pS, struct easyGene *ret)
/* Create a easyGene out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new easyGene */
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
ret->feat = sqlStringComma(&s);
ret->R = sqlFloatComma(&s);
ret->frame = sqlUnsignedComma(&s);
ret->orf = sqlStringComma(&s);
ret->startCodon = sqlStringComma(&s);
ret->logOdds = sqlFloatComma(&s);
ret->descriptor = sqlStringComma(&s);
ret->swissProt = sqlStringComma(&s);
sqlFixedStringComma(&s, ret->genbank, sizeof(ret->genbank));
*pS = s;
return ret;
}

void easyGeneFree(struct easyGene **pEl)
/* Free a single dynamically allocated easyGene such as created
 * with easyGeneLoad(). */
{
struct easyGene *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->feat);
freeMem(el->orf);
freeMem(el->startCodon);
freeMem(el->descriptor);
freeMem(el->swissProt);
freez(pEl);
}

void easyGeneFreeList(struct easyGene **pList)
/* Free a list of dynamically allocated easyGene's */
{
struct easyGene *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    easyGeneFree(&el);
    }
*pList = NULL;
}

void easyGeneOutput(struct easyGene *el, FILE *f, char sep, char lastSep) 
/* Print out easyGene.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%s", el->feat);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%g", el->R);
fputc(sep,f);
fprintf(f, "%u", el->frame);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->orf);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->startCodon);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%g", el->logOdds);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->descriptor);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->swissProt);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->genbank);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

