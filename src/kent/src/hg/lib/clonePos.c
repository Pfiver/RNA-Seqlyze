/* clonePos.c was originally generated by the autoSql program, which also 
 * generated clonePos.h and clonePos.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "clonePos.h"


void clonePosStaticLoad(char **row, struct clonePos *ret)
/* Load a row from clonePos table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->name = row[0];
ret->seqSize = sqlUnsigned(row[1]);
ret->phase = sqlUnsigned(row[2]);
ret->chrom = row[3];
ret->chromStart = sqlUnsigned(row[4]);
ret->chromEnd = sqlUnsigned(row[5]);
safecpy(ret->stage, sizeof(ret->stage), row[6]);
ret->faFile = row[7];
}

struct clonePos *clonePosLoad(char **row)
/* Load a clonePos from row fetched with select * from clonePos
 * from database.  Dispose of this with clonePosFree(). */
{
struct clonePos *ret;

AllocVar(ret);
ret->name = cloneString(row[0]);
ret->seqSize = sqlUnsigned(row[1]);
ret->phase = sqlUnsigned(row[2]);
ret->chrom = cloneString(row[3]);
ret->chromStart = sqlUnsigned(row[4]);
ret->chromEnd = sqlUnsigned(row[5]);
safecpy(ret->stage, sizeof(ret->stage), row[6]);
ret->faFile = cloneString(row[7]);
return ret;
}

struct clonePos *clonePosLoadAll(char *fileName) 
/* Load all clonePos from a whitespace-separated file.
 * Dispose of this with clonePosFreeList(). */
{
struct clonePos *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[8];

while (lineFileRow(lf, row))
    {
    el = clonePosLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct clonePos *clonePosLoadAllByChar(char *fileName, char chopper) 
/* Load all clonePos from a chopper separated file.
 * Dispose of this with clonePosFreeList(). */
{
struct clonePos *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[8];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = clonePosLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct clonePos *clonePosCommaIn(char **pS, struct clonePos *ret)
/* Create a clonePos out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new clonePos */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->seqSize = sqlUnsignedComma(&s);
ret->phase = sqlUnsignedComma(&s);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->stage, sizeof(ret->stage));
ret->faFile = sqlStringComma(&s);
*pS = s;
return ret;
}

void clonePosFree(struct clonePos **pEl)
/* Free a single dynamically allocated clonePos such as created
 * with clonePosLoad(). */
{
struct clonePos *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->chrom);
freeMem(el->faFile);
freez(pEl);
}

void clonePosFreeList(struct clonePos **pList)
/* Free a list of dynamically allocated clonePos's */
{
struct clonePos *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    clonePosFree(&el);
    }
*pList = NULL;
}

void clonePosOutput(struct clonePos *el, FILE *f, char sep, char lastSep) 
/* Print out clonePos.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->seqSize);
fputc(sep,f);
fprintf(f, "%u", el->phase);
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
fprintf(f, "%s", el->stage);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->faFile);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

