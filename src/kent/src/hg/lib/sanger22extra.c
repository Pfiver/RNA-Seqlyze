/* sanger22extra.c was originally generated by the autoSql program, which also 
 * generated sanger22extra.h and sanger22extra.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "jksql.h"
#include "sanger22extra.h"


void sanger22extraStaticLoad(char **row, struct sanger22extra *ret)
/* Load a row from sanger22extra table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->name = row[0];
ret->locus = row[1];
ret->description = row[2];
ret->geneType = row[3];
ret->cdsType = row[4];
}

struct sanger22extra *sanger22extraLoad(char **row)
/* Load a sanger22extra from row fetched with select * from sanger22extra
 * from database.  Dispose of this with sanger22extraFree(). */
{
struct sanger22extra *ret;

AllocVar(ret);
ret->name = cloneString(row[0]);
ret->locus = cloneString(row[1]);
ret->description = cloneString(row[2]);
ret->geneType = cloneString(row[3]);
ret->cdsType = cloneString(row[4]);
return ret;
}

struct sanger22extra *sanger22extraLoadAll(char *fileName) 
/* Load all sanger22extra from a tab-separated file.
 * Dispose of this with sanger22extraFreeList(). */
{
struct sanger22extra *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[5];

while (lineFileRow(lf, row))
    {
    el = sanger22extraLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct sanger22extra *sanger22extraCommaIn(char **pS, struct sanger22extra *ret)
/* Create a sanger22extra out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new sanger22extra */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->locus = sqlStringComma(&s);
ret->description = sqlStringComma(&s);
ret->geneType = sqlStringComma(&s);
ret->cdsType = sqlStringComma(&s);
*pS = s;
return ret;
}

void sanger22extraFree(struct sanger22extra **pEl)
/* Free a single dynamically allocated sanger22extra such as created
 * with sanger22extraLoad(). */
{
struct sanger22extra *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->locus);
freeMem(el->description);
freeMem(el->geneType);
freeMem(el->cdsType);
freez(pEl);
}

void sanger22extraFreeList(struct sanger22extra **pList)
/* Free a list of dynamically allocated sanger22extra's */
{
struct sanger22extra *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    sanger22extraFree(&el);
    }
*pList = NULL;
}

void sanger22extraOutput(struct sanger22extra *el, FILE *f, char sep, char lastSep) 
/* Print out sanger22extra.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->locus);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->description);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->geneType);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->cdsType);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* ---------------- End of AutoSQL generated code. ----------------------- */

char *sanger22extraCreate = 
"CREATE TABLE %s (\n"
"    name varchar(255) not null,	# Transcript name\n"
"    locus varchar(255) not null,	# Possibly biological short name\n"
"    description longblob not null,	# Description from Sanger gene GFFs\n"
"    geneType varchar(255) not null,	# Type field from Sanger gene GFFs\n"
"    cdsType varchar(255) not null,	# Type field from Sanger CDS GFFs\n"
"              #Indices\n"
"    PRIMARY KEY(name)\n"
")\n";

