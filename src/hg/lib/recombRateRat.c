/* recombRateRat.c was originally generated by the autoSql program, which also 
 * generated recombRateRat.h and recombRateRat.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "recombRateRat.h"


void recombRateRatStaticLoad(char **row, struct recombRateRat *ret)
/* Load a row from recombRateRat table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->shrspAvg = atof(row[4]);
ret->fhhAvg = atof(row[5]);
}

struct recombRateRat *recombRateRatLoad(char **row)
/* Load a recombRateRat from row fetched with select * from recombRateRat
 * from database.  Dispose of this with recombRateRatFree(). */
{
struct recombRateRat *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->shrspAvg = atof(row[4]);
ret->fhhAvg = atof(row[5]);
return ret;
}

struct recombRateRat *recombRateRatLoadAll(char *fileName) 
/* Load all recombRateRat from a whitespace-separated file.
 * Dispose of this with recombRateRatFreeList(). */
{
struct recombRateRat *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[6];

while (lineFileRow(lf, row))
    {
    el = recombRateRatLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct recombRateRat *recombRateRatLoadAllByChar(char *fileName, char chopper) 
/* Load all recombRateRat from a chopper separated file.
 * Dispose of this with recombRateRatFreeList(). */
{
struct recombRateRat *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[6];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = recombRateRatLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct recombRateRat *recombRateRatCommaIn(char **pS, struct recombRateRat *ret)
/* Create a recombRateRat out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new recombRateRat */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->shrspAvg = sqlFloatComma(&s);
ret->fhhAvg = sqlFloatComma(&s);
*pS = s;
return ret;
}

void recombRateRatFree(struct recombRateRat **pEl)
/* Free a single dynamically allocated recombRateRat such as created
 * with recombRateRatLoad(). */
{
struct recombRateRat *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void recombRateRatFreeList(struct recombRateRat **pList)
/* Free a list of dynamically allocated recombRateRat's */
{
struct recombRateRat *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    recombRateRatFree(&el);
    }
*pList = NULL;
}

void recombRateRatOutput(struct recombRateRat *el, FILE *f, char sep, char lastSep) 
/* Print out recombRateRat.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%f", el->shrspAvg);
fputc(sep,f);
fprintf(f, "%f", el->fhhAvg);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

