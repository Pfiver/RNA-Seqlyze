/* cdsEvidence.c was originally generated by the autoSql program, which also 
 * generated cdsEvidence.h and cdsEvidence.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "bed.h"
#include "cdsEvidence.h"


struct cdsEvidence *cdsEvidenceLoad(char **row)
/* Load a cdsEvidence from row fetched with select * from cdsEvidence
 * from database.  Dispose of this with cdsEvidenceFree(). */
{
struct cdsEvidence *ret;

AllocVar(ret);
ret->cdsCount = sqlSigned(row[8]);
ret->name = cloneString(row[0]);
ret->start = sqlSigned(row[1]);
ret->end = sqlSigned(row[2]);
ret->source = cloneString(row[3]);
ret->accession = cloneString(row[4]);
ret->score = sqlDouble(row[5]);
ret->startComplete = sqlUnsigned(row[6]);
ret->endComplete = sqlUnsigned(row[7]);
{
int sizeOne;
sqlSignedDynamicArray(row[9], &ret->cdsStarts, &sizeOne);
assert(sizeOne == ret->cdsCount);
}
{
int sizeOne;
sqlSignedDynamicArray(row[10], &ret->cdsSizes, &sizeOne);
assert(sizeOne == ret->cdsCount);
}
return ret;
}

struct cdsEvidence *cdsEvidenceLoadAll(char *fileName) 
/* Load all cdsEvidence from a whitespace-separated file.
 * Dispose of this with cdsEvidenceFreeList(). */
{
struct cdsEvidence *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[11];

while (lineFileRow(lf, row))
    {
    el = cdsEvidenceLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct cdsEvidence *cdsEvidenceLoadAllByChar(char *fileName, char chopper) 
/* Load all cdsEvidence from a chopper separated file.
 * Dispose of this with cdsEvidenceFreeList(). */
{
struct cdsEvidence *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[11];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = cdsEvidenceLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct cdsEvidence *cdsEvidenceCommaIn(char **pS, struct cdsEvidence *ret)
/* Create a cdsEvidence out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdsEvidence */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->start = sqlSignedComma(&s);
ret->end = sqlSignedComma(&s);
ret->source = sqlStringComma(&s);
ret->accession = sqlStringComma(&s);
ret->score = sqlDoubleComma(&s);
ret->startComplete = sqlUnsignedComma(&s);
ret->endComplete = sqlUnsignedComma(&s);
ret->cdsCount = sqlSignedComma(&s);
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->cdsStarts, ret->cdsCount);
for (i=0; i<ret->cdsCount; ++i)
    {
    ret->cdsStarts[i] = sqlSignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->cdsSizes, ret->cdsCount);
for (i=0; i<ret->cdsCount; ++i)
    {
    ret->cdsSizes[i] = sqlSignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
*pS = s;
return ret;
}

void cdsEvidenceFree(struct cdsEvidence **pEl)
/* Free a single dynamically allocated cdsEvidence such as created
 * with cdsEvidenceLoad(). */
{
struct cdsEvidence *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->source);
freeMem(el->accession);
freeMem(el->cdsStarts);
freeMem(el->cdsSizes);
freez(pEl);
}

void cdsEvidenceFreeList(struct cdsEvidence **pList)
/* Free a list of dynamically allocated cdsEvidence's */
{
struct cdsEvidence *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    cdsEvidenceFree(&el);
    }
*pList = NULL;
}

void cdsEvidenceOutput(struct cdsEvidence *el, FILE *f, char sep, char lastSep) 
/* Print out cdsEvidence.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->start);
fputc(sep,f);
fprintf(f, "%d", el->end);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->source);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->accession);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%g", el->score);
fputc(sep,f);
fprintf(f, "%u", el->startComplete);
fputc(sep,f);
fprintf(f, "%u", el->endComplete);
fputc(sep,f);
fprintf(f, "%d", el->cdsCount);
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->cdsCount; ++i)
    {
    fprintf(f, "%d", el->cdsStarts[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->cdsCount; ++i)
    {
    fprintf(f, "%d", el->cdsSizes[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

struct hash *cdsEvidenceReadAllIntoHash(char *fileName)
/* Return hash full of cdsEvidence keyed by transcript name. */
{
struct lineFile *lf = lineFileOpen(fileName, TRUE);
struct hash *hash = hashNew(18);
char *row[CDSEVIDENCE_NUM_COLS];
while (lineFileRowTab(lf, row))
    {
    struct cdsEvidence *cds = cdsEvidenceLoad(row);
    if (hashLookup(hash, cds->name))
        errAbort("%s duplicated in %s, perhaps you want to run txCdsPick?",
		cds->name, fileName);
    hashAdd(hash, cds->name, cds);
    }
lineFileClose(&lf);
return hash;
}

void cdsEvidenceSetBedThick(struct cdsEvidence *cds, struct bed *bed,
			    const boolean freeOfCdsErrors)
/* Set thickStart/thickEnd on bed from cdsEvidence. */
{
if (cds == NULL || !freeOfCdsErrors)
    {
    bed->thickStart = bed->thickEnd = bed->chromStart;
    return;
    }
int txCdsStart = cds->start, txCdsEnd = cds->end;
if (bed->strand[0] == '-')
    {
    int txSize = bedTotalBlockSize(bed);
    reverseIntRange(&txCdsStart, &txCdsEnd, txSize);
    }
int i;
int txStart = 0, txEnd;
for (i=0; i<bed->blockCount; ++i)
    {
    int blockSize = bed->blockSizes[i];
    int exonStart = bed->chromStarts[i] + bed->chromStart;
    txEnd = txStart + blockSize;
    if (txStart <= txCdsStart && txCdsStart < txEnd)
        {
	int offset = txCdsStart - txStart;
	bed->thickStart = exonStart + offset;
	}
    if (txStart < txCdsEnd && txCdsEnd <= txEnd)
        {
	int offset = txCdsEnd - txStart;
	bed->thickEnd = exonStart + offset;
	}
    txStart = txEnd;
    }
}

int cdsEvidenceCmpScore(const void *va, const void *vb)
/* Compare to sort based on score (descending). */
{
const struct cdsEvidence *a = *((struct cdsEvidence **)va);
const struct cdsEvidence *b = *((struct cdsEvidence **)vb);
double diff = b->score - a->score;
if (diff < 0)
    return -1;
else if (diff > 0)
    return 1;
else
    return 0;
}

