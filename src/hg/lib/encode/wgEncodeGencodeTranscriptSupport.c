/* wgEncodeGencodeTranscriptSupport.c was originally generated by the autoSql program, which also 
 * generated wgEncodeGencodeTranscriptSupport.h and wgEncodeGencodeTranscriptSupport.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "encode/wgEncodeGencodeTranscriptSupport.h"


void wgEncodeGencodeTranscriptSupportStaticLoad(char **row, struct wgEncodeGencodeTranscriptSupport *ret)
/* Load a row from wgEncodeGencodeTranscriptSupport table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->transcriptId = row[0];
ret->seqId = row[1];
ret->seqSrc = row[2];
}

struct wgEncodeGencodeTranscriptSupport *wgEncodeGencodeTranscriptSupportLoad(char **row)
/* Load a wgEncodeGencodeTranscriptSupport from row fetched with select * from wgEncodeGencodeTranscriptSupport
 * from database.  Dispose of this with wgEncodeGencodeTranscriptSupportFree(). */
{
struct wgEncodeGencodeTranscriptSupport *ret;

AllocVar(ret);
ret->transcriptId = cloneString(row[0]);
ret->seqId = cloneString(row[1]);
ret->seqSrc = cloneString(row[2]);
return ret;
}

struct wgEncodeGencodeTranscriptSupport *wgEncodeGencodeTranscriptSupportLoadAll(char *fileName) 
/* Load all wgEncodeGencodeTranscriptSupport from a whitespace-separated file.
 * Dispose of this with wgEncodeGencodeTranscriptSupportFreeList(). */
{
struct wgEncodeGencodeTranscriptSupport *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileRow(lf, row))
    {
    el = wgEncodeGencodeTranscriptSupportLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeTranscriptSupport *wgEncodeGencodeTranscriptSupportLoadAllByChar(char *fileName, char chopper) 
/* Load all wgEncodeGencodeTranscriptSupport from a chopper separated file.
 * Dispose of this with wgEncodeGencodeTranscriptSupportFreeList(). */
{
struct wgEncodeGencodeTranscriptSupport *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = wgEncodeGencodeTranscriptSupportLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeTranscriptSupport *wgEncodeGencodeTranscriptSupportCommaIn(char **pS, struct wgEncodeGencodeTranscriptSupport *ret)
/* Create a wgEncodeGencodeTranscriptSupport out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wgEncodeGencodeTranscriptSupport */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->transcriptId = sqlStringComma(&s);
ret->seqId = sqlStringComma(&s);
ret->seqSrc = sqlStringComma(&s);
*pS = s;
return ret;
}

void wgEncodeGencodeTranscriptSupportFree(struct wgEncodeGencodeTranscriptSupport **pEl)
/* Free a single dynamically allocated wgEncodeGencodeTranscriptSupport such as created
 * with wgEncodeGencodeTranscriptSupportLoad(). */
{
struct wgEncodeGencodeTranscriptSupport *el;

if ((el = *pEl) == NULL) return;
freeMem(el->transcriptId);
freeMem(el->seqId);
freeMem(el->seqSrc);
freez(pEl);
}

void wgEncodeGencodeTranscriptSupportFreeList(struct wgEncodeGencodeTranscriptSupport **pList)
/* Free a list of dynamically allocated wgEncodeGencodeTranscriptSupport's */
{
struct wgEncodeGencodeTranscriptSupport *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    wgEncodeGencodeTranscriptSupportFree(&el);
    }
*pList = NULL;
}

void wgEncodeGencodeTranscriptSupportOutput(struct wgEncodeGencodeTranscriptSupport *el, FILE *f, char sep, char lastSep) 
/* Print out wgEncodeGencodeTranscriptSupport.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->transcriptId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->seqId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->seqSrc);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

