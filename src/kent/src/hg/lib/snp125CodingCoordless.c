/* snp125CodingCoordless.c was originally generated by the autoSql program, which also 
 * generated snp125CodingCoordless.h and snp125CodingCoordless.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "snp125CodingCoordless.h"


/* definitions for frame column */
static char *values_frame[] = {"1", "2", "3", NULL};
static struct hash *valhash_frame = NULL;

struct snp125CodingCoordless *snp125CodingCoordlessLoad(char **row)
/* Load a snp125CodingCoordless from row fetched with select * from snp125CodingCoordless
 * from database.  Dispose of this with snp125CodingCoordlessFree(). */
{
struct snp125CodingCoordless *ret;

AllocVar(ret);
ret->alleleCount = sqlSigned(row[3]);
ret->name = cloneString(row[0]);
ret->transcript = cloneString(row[1]);
ret->frame = sqlEnumParse(row[2], values_frame, &valhash_frame);
{
int sizeOne;
sqlUshortDynamicArray(row[4], &ret->funcCodes, &sizeOne);
assert(sizeOne == ret->alleleCount);
}
{
int sizeOne;
sqlStringDynamicArray(row[5], &ret->alleles, &sizeOne);
assert(sizeOne == ret->alleleCount);
}
{
int sizeOne;
sqlStringDynamicArray(row[6], &ret->codons, &sizeOne);
assert(sizeOne == ret->alleleCount);
}
{
int sizeOne;
sqlStringDynamicArray(row[7], &ret->peptides, &sizeOne);
assert(sizeOne == ret->alleleCount);
}
return ret;
}

struct snp125CodingCoordless *snp125CodingCoordlessLoadAll(char *fileName) 
/* Load all snp125CodingCoordless from a whitespace-separated file.
 * Dispose of this with snp125CodingCoordlessFreeList(). */
{
struct snp125CodingCoordless *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[8];

while (lineFileRow(lf, row))
    {
    el = snp125CodingCoordlessLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct snp125CodingCoordless *snp125CodingCoordlessLoadAllByChar(char *fileName, char chopper) 
/* Load all snp125CodingCoordless from a chopper separated file.
 * Dispose of this with snp125CodingCoordlessFreeList(). */
{
struct snp125CodingCoordless *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[8];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = snp125CodingCoordlessLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct snp125CodingCoordless *snp125CodingCoordlessCommaIn(char **pS, struct snp125CodingCoordless *ret)
/* Create a snp125CodingCoordless out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new snp125CodingCoordless */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->transcript = sqlStringComma(&s);
ret->frame = sqlEnumComma(&s, values_frame, &valhash_frame);
ret->alleleCount = sqlSignedComma(&s);
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->funcCodes, ret->alleleCount);
for (i=0; i<ret->alleleCount; ++i)
    {
    ret->funcCodes[i] = sqlUnsignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->alleles, ret->alleleCount);
for (i=0; i<ret->alleleCount; ++i)
    {
    ret->alleles[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->codons, ret->alleleCount);
for (i=0; i<ret->alleleCount; ++i)
    {
    ret->codons[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->peptides, ret->alleleCount);
for (i=0; i<ret->alleleCount; ++i)
    {
    ret->peptides[i] = sqlStringComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
*pS = s;
return ret;
}

void snp125CodingCoordlessFree(struct snp125CodingCoordless **pEl)
/* Free a single dynamically allocated snp125CodingCoordless such as created
 * with snp125CodingCoordlessLoad(). */
{
struct snp125CodingCoordless *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->transcript);
freeMem(el->funcCodes);
/* All strings in alleles are allocated at once, so only need to free first. */
if (el->alleles != NULL)
    freeMem(el->alleles[0]);
freeMem(el->alleles);
/* All strings in codons are allocated at once, so only need to free first. */
if (el->codons != NULL)
    freeMem(el->codons[0]);
freeMem(el->codons);
/* All strings in peptides are allocated at once, so only need to free first. */
if (el->peptides != NULL)
    freeMem(el->peptides[0]);
freeMem(el->peptides);
freez(pEl);
}

void snp125CodingCoordlessFreeList(struct snp125CodingCoordless **pList)
/* Free a list of dynamically allocated snp125CodingCoordless's */
{
struct snp125CodingCoordless *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    snp125CodingCoordlessFree(&el);
    }
*pList = NULL;
}

void snp125CodingCoordlessOutput(struct snp125CodingCoordless *el, FILE *f, char sep, char lastSep) 
/* Print out snp125CodingCoordless.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->transcript);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
sqlEnumPrint(f, el->frame, values_frame);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->alleleCount);
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->alleleCount; ++i)
    {
    fprintf(f, "%u", el->funcCodes[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->alleleCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->alleles[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->alleleCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->codons[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->alleleCount; ++i)
    {
    if (sep == ',') fputc('"',f);
    fprintf(f, "%s", el->peptides[i]);
    if (sep == ',') fputc('"',f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

