/* delConrad2.c was originally generated by the autoSql program, which also 
 * generated delConrad2.h and delConrad2.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "delConrad2.h"


void delConrad2StaticLoad(char **row, struct delConrad2 *ret)
/* Load a row from delConrad2 table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
ret->count1 = sqlUnsigned(row[8]);
ret->count2 = sqlUnsigned(row[9]);
ret->offspring = row[10];
ret->population = row[11];
}

struct delConrad2 *delConrad2Load(char **row)
/* Load a delConrad2 from row fetched with select * from delConrad2
 * from database.  Dispose of this with delConrad2Free(). */
{
struct delConrad2 *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
ret->count1 = sqlUnsigned(row[8]);
ret->count2 = sqlUnsigned(row[9]);
ret->offspring = cloneString(row[10]);
ret->population = cloneString(row[11]);
return ret;
}

struct delConrad2 *delConrad2LoadAll(char *fileName) 
/* Load all delConrad2 from a whitespace-separated file.
 * Dispose of this with delConrad2FreeList(). */
{
struct delConrad2 *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[12];

while (lineFileRow(lf, row))
    {
    el = delConrad2Load(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct delConrad2 *delConrad2LoadAllByChar(char *fileName, char chopper) 
/* Load all delConrad2 from a chopper separated file.
 * Dispose of this with delConrad2FreeList(). */
{
struct delConrad2 *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[12];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = delConrad2Load(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct delConrad2 *delConrad2CommaIn(char **pS, struct delConrad2 *ret)
/* Create a delConrad2 out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new delConrad2 */
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
ret->thickStart = sqlUnsignedComma(&s);
ret->thickEnd = sqlUnsignedComma(&s);
ret->count1 = sqlUnsignedComma(&s);
ret->count2 = sqlUnsignedComma(&s);
ret->offspring = sqlStringComma(&s);
ret->population = sqlStringComma(&s);
*pS = s;
return ret;
}

void delConrad2Free(struct delConrad2 **pEl)
/* Free a single dynamically allocated delConrad2 such as created
 * with delConrad2Load(). */
{
struct delConrad2 *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->offspring);
freeMem(el->population);
freez(pEl);
}

void delConrad2FreeList(struct delConrad2 **pList)
/* Free a list of dynamically allocated delConrad2's */
{
struct delConrad2 *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    delConrad2Free(&el);
    }
*pList = NULL;
}

void delConrad2Output(struct delConrad2 *el, FILE *f, char sep, char lastSep) 
/* Print out delConrad2.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%u", el->thickStart);
fputc(sep,f);
fprintf(f, "%u", el->thickEnd);
fputc(sep,f);
fprintf(f, "%u", el->count1);
fputc(sep,f);
fprintf(f, "%u", el->count2);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->offspring);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->population);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

