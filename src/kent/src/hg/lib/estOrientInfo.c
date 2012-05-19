/* estOrientInfo.c was originally generated by the autoSql program, which also 
 * generated estOrientInfo.h and estOrientInfo.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "localmem.h"
#include "estOrientInfo.h"


static char *createString = 
    "CREATE TABLE %s (\n"
    "    bin smallint unsigned not null,        # Bin for fast index\n"
    "    chrom varchar(255) not null,	        # Human chromosome or FPC contig\n"
    "    chromStart int unsigned not null,	# Start position in chromosome\n"
    "    chromEnd int unsigned not null,	# End position in chromosome\n"
    "    name varchar(255) not null,	        # Accession of EST\n"
    "    intronOrientation smallint not null,	# Orientation of introns with respect to EST\n"
    "    sizePolyA smallint not null,	        # Number of trailing A's\n"
    "    revSizePolyA smallint not null,	# Number of trailing A's on reverse strand\n"
    "    signalPos smallint not null,	        # Position of start of polyA signal relative to end of EST or 0 if no signal\n"
    "    revSignalPos smallint not null,	# PolyA signal position on reverse strand if any\n"
    "              #Indices\n"
    "    INDEX(chrom(%d),bin),\n"
    "    INDEX(chrom(%d),chromStart),\n"
    "    INDEX(chrom(%d),chromEnd),\n"
    "    INDEX(name(20))\n"
    ")\n";

void estOrientInfoStaticLoad(char **row, struct estOrientInfo *ret)
/* Load a row from estOrientInfo table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->intronOrientation = sqlSigned(row[4]);
ret->sizePolyA = sqlSigned(row[5]);
ret->revSizePolyA = sqlSigned(row[6]);
ret->signalPos = sqlSigned(row[7]);
ret->revSignalPos = sqlSigned(row[8]);
}

struct estOrientInfo *estOrientInfoLoad(char **row)
/* Load a estOrientInfo from row fetched with select * from estOrientInfo
 * from database.  Dispose of this with estOrientInfoFree(). */
{
struct estOrientInfo *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->intronOrientation = sqlSigned(row[4]);
ret->sizePolyA = sqlSigned(row[5]);
ret->revSizePolyA = sqlSigned(row[6]);
ret->signalPos = sqlSigned(row[7]);
ret->revSignalPos = sqlSigned(row[8]);
return ret;
}

struct estOrientInfo *estOrientInfoLoadLm(char **row, struct lm *lm)
/* Load a estOrientInfo row into local memory struct. */
{
struct estOrientInfo *ret;

lmAllocVar(lm, ret);
ret->chrom = lmCloneString(lm, row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->intronOrientation = sqlSigned(row[4]);
ret->sizePolyA = sqlSigned(row[5]);
ret->revSizePolyA = sqlSigned(row[6]);
ret->signalPos = sqlSigned(row[7]);
ret->revSignalPos = sqlSigned(row[8]);
return ret;
}

struct estOrientInfo *estOrientInfoLoadAll(char *fileName) 
/* Load all estOrientInfo from a tab-separated file.
 * Dispose of this with estOrientInfoFreeList(). */
{
struct estOrientInfo *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[9];

while (lineFileRow(lf, row))
    {
    el = estOrientInfoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct estOrientInfo *estOrientInfoLoadWhere(struct sqlConnection *conn, char *table, char *where)
/* Load all estOrientInfo from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with estOrientInfoFreeList(). */
{
struct estOrientInfo *list = NULL, *el;
struct dyString *query = dyStringNew(256);
struct sqlResult *sr;
char **row;

dyStringPrintf(query, "select * from %s", table);
if (where != NULL)
    dyStringPrintf(query, " where %s", where);
sr = sqlGetResult(conn, query->string);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = estOrientInfoLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
dyStringFree(&query);
return list;
}

struct estOrientInfo *estOrientInfoCommaIn(char **pS, struct estOrientInfo *ret)
/* Create a estOrientInfo out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new estOrientInfo */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->intronOrientation = sqlSignedComma(&s);
ret->sizePolyA = sqlSignedComma(&s);
ret->revSizePolyA = sqlSignedComma(&s);
ret->signalPos = sqlSignedComma(&s);
ret->revSignalPos = sqlSignedComma(&s);
*pS = s;
return ret;
}

void estOrientInfoFree(struct estOrientInfo **pEl)
/* Free a single dynamically allocated estOrientInfo such as created
 * with estOrientInfoLoad(). */
{
struct estOrientInfo *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void estOrientInfoFreeList(struct estOrientInfo **pList)
/* Free a list of dynamically allocated estOrientInfo's */
{
struct estOrientInfo *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    estOrientInfoFree(&el);
    }
*pList = NULL;
}

void estOrientInfoOutput(struct estOrientInfo *el, FILE *f, char sep, char lastSep) 
/* Print out estOrientInfo.  Separate fields with sep. Follow last field with lastSep. */
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
fprintf(f, "%d", el->intronOrientation);
fputc(sep,f);
fprintf(f, "%d", el->sizePolyA);
fputc(sep,f);
fprintf(f, "%d", el->revSizePolyA);
fputc(sep,f);
fprintf(f, "%d", el->signalPos);
fputc(sep,f);
fprintf(f, "%d", el->revSignalPos);
fputc(lastSep,f);
}

char *estOrientInfoGetCreateSql(char *table, int chromIdxLen)
/* Get SQL to create an estOrientInfo table. chromIdxLen is the number of
 * chars at that start of chrom to use for the index. */
{
struct dyString *sqlCmd = newDyString(2048);
char *sqlCmdStr;
dyStringPrintf(sqlCmd, createString, table, chromIdxLen, chromIdxLen, chromIdxLen);
sqlCmdStr = cloneString(sqlCmd->string);
dyStringFree(&sqlCmd);
return sqlCmdStr;
}

