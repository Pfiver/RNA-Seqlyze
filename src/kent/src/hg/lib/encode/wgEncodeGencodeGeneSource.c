/* wgEncodeGencodeGeneSource.c was originally generated by the autoSql program, which also 
 * generated wgEncodeGencodeGeneSource.h and wgEncodeGencodeGeneSource.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "encode/wgEncodeGencodeGeneSource.h"


void wgEncodeGencodeGeneSourceStaticLoad(char **row, struct wgEncodeGencodeGeneSource *ret)
/* Load a row from wgEncodeGencodeGeneSource table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->geneId = row[0];
ret->source = row[1];
}

struct wgEncodeGencodeGeneSource *wgEncodeGencodeGeneSourceLoad(char **row)
/* Load a wgEncodeGencodeGeneSource from row fetched with select * from wgEncodeGencodeGeneSource
 * from database.  Dispose of this with wgEncodeGencodeGeneSourceFree(). */
{
struct wgEncodeGencodeGeneSource *ret;

AllocVar(ret);
ret->geneId = cloneString(row[0]);
ret->source = cloneString(row[1]);
return ret;
}

struct wgEncodeGencodeGeneSource *wgEncodeGencodeGeneSourceLoadAll(char *fileName) 
/* Load all wgEncodeGencodeGeneSource from a whitespace-separated file.
 * Dispose of this with wgEncodeGencodeGeneSourceFreeList(). */
{
struct wgEncodeGencodeGeneSource *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileRow(lf, row))
    {
    el = wgEncodeGencodeGeneSourceLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeGeneSource *wgEncodeGencodeGeneSourceLoadAllByChar(char *fileName, char chopper) 
/* Load all wgEncodeGencodeGeneSource from a chopper separated file.
 * Dispose of this with wgEncodeGencodeGeneSourceFreeList(). */
{
struct wgEncodeGencodeGeneSource *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[2];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = wgEncodeGencodeGeneSourceLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeGeneSource *wgEncodeGencodeGeneSourceCommaIn(char **pS, struct wgEncodeGencodeGeneSource *ret)
/* Create a wgEncodeGencodeGeneSource out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wgEncodeGencodeGeneSource */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->geneId = sqlStringComma(&s);
ret->source = sqlStringComma(&s);
*pS = s;
return ret;
}

void wgEncodeGencodeGeneSourceFree(struct wgEncodeGencodeGeneSource **pEl)
/* Free a single dynamically allocated wgEncodeGencodeGeneSource such as created
 * with wgEncodeGencodeGeneSourceLoad(). */
{
struct wgEncodeGencodeGeneSource *el;

if ((el = *pEl) == NULL) return;
freeMem(el->geneId);
freeMem(el->source);
freez(pEl);
}

void wgEncodeGencodeGeneSourceFreeList(struct wgEncodeGencodeGeneSource **pList)
/* Free a list of dynamically allocated wgEncodeGencodeGeneSource's */
{
struct wgEncodeGencodeGeneSource *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    wgEncodeGencodeGeneSourceFree(&el);
    }
*pList = NULL;
}

void wgEncodeGencodeGeneSourceOutput(struct wgEncodeGencodeGeneSource *el, FILE *f, char sep, char lastSep) 
/* Print out wgEncodeGencodeGeneSource.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->geneId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->source);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

