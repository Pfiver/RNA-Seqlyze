/* targetDb.c was originally generated by the autoSql program, which also 
 * generated targetDb.h and targetDb.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "targetDb.h"


void targetDbStaticLoad(char **row, struct targetDb *ret)
/* Load a row from targetDb table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->name = row[0];
ret->description = row[1];
ret->db = row[2];
ret->pslTable = row[3];
ret->seqTable = row[4];
ret->extFileTable = row[5];
ret->seqFile = row[6];
ret->priority = sqlFloat(row[7]);
ret->time = row[8];
ret->settings = row[9];
}

struct targetDb *targetDbLoad(char **row)
/* Load a targetDb from row fetched with select * from targetDb
 * from database.  Dispose of this with targetDbFree(). */
{
struct targetDb *ret;

AllocVar(ret);
ret->name = cloneString(row[0]);
ret->description = cloneString(row[1]);
ret->db = cloneString(row[2]);
ret->pslTable = cloneString(row[3]);
ret->seqTable = cloneString(row[4]);
ret->extFileTable = cloneString(row[5]);
ret->seqFile = cloneString(row[6]);
ret->priority = sqlFloat(row[7]);
ret->time = cloneString(row[8]);
ret->settings = cloneString(row[9]);
return ret;
}

struct targetDb *targetDbLoadAll(char *fileName) 
/* Load all targetDb from a whitespace-separated file.
 * Dispose of this with targetDbFreeList(). */
{
struct targetDb *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[10];

while (lineFileRow(lf, row))
    {
    el = targetDbLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct targetDb *targetDbLoadAllByChar(char *fileName, char chopper) 
/* Load all targetDb from a chopper separated file.
 * Dispose of this with targetDbFreeList(). */
{
struct targetDb *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[10];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = targetDbLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct targetDb *targetDbCommaIn(char **pS, struct targetDb *ret)
/* Create a targetDb out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new targetDb */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->description = sqlStringComma(&s);
ret->db = sqlStringComma(&s);
ret->pslTable = sqlStringComma(&s);
ret->seqTable = sqlStringComma(&s);
ret->extFileTable = sqlStringComma(&s);
ret->seqFile = sqlStringComma(&s);
ret->priority = sqlFloatComma(&s);
ret->time = sqlStringComma(&s);
ret->settings = sqlStringComma(&s);
*pS = s;
return ret;
}

void targetDbFree(struct targetDb **pEl)
/* Free a single dynamically allocated targetDb such as created
 * with targetDbLoad(). */
{
struct targetDb *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->description);
freeMem(el->db);
freeMem(el->pslTable);
freeMem(el->seqTable);
freeMem(el->extFileTable);
freeMem(el->seqFile);
freeMem(el->time);
freeMem(el->settings);
freez(pEl);
}

void targetDbFreeList(struct targetDb **pList)
/* Free a list of dynamically allocated targetDb's */
{
struct targetDb *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    targetDbFree(&el);
    }
*pList = NULL;
}

void targetDbOutput(struct targetDb *el, FILE *f, char sep, char lastSep) 
/* Print out targetDb.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->description);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->db);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->pslTable);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->seqTable);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->extFileTable);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->seqFile);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%g", el->priority);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->time);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->settings);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#include "portable.h"
#include "hdb.h"
#include "ra.h"

static boolean timeMoreRecentThanTable(int time, struct sqlConnection *conn,
				char *table)
/* Return TRUE if the given UNIX time is more recent than the time that 
 * table was last updated. */
{
if (! sqlTableExists(conn, table))
    return FALSE;
int tableUpdateTime = sqlTableUpdateTime(conn, table);
return (time > tableUpdateTime);
}

static boolean timeMoreRecentThanFile(int time, char *fileName)
/* Return TRUE if the given UNIX time is more recent than the time that
 * fileName was last modified. */
{
if (! fileExists(fileName))
    return FALSE;
int fileUpdateTime = fileModTime(fileName);
return (time > fileUpdateTime);
}

struct targetDb *targetDbMaybeLoad(struct sqlConnection *conn, char **row)
/* If row specifies a target whose tables and file exist, and are not newer
 * than target, allocate and return a targetDb; otherwise, return NULL
 * and log a warning to stderr for QA monitoring. */
{
struct targetDb target;
targetDbStaticLoad(row, &target);
int time = sqlDateToUnixTime(target.time);
if (timeMoreRecentThanTable(time, conn, target.pslTable) &&
    (isEmpty(target.seqTable) ||
     timeMoreRecentThanTable(time, conn, target.seqTable)) &&
    (isEmpty(target.extFileTable) ||
     timeMoreRecentThanTable(time, conn, target.extFileTable)) &&
/* Here is where we could add a configurable slush factor: */
    timeMoreRecentThanFile(time, target.seqFile))
    {
    return targetDbLoad(row);
    }
else
    fprintf(stderr, "targetDb entry %s is dated %s -- older than at "
	    "least one of its db tables (%s, %s, %s) "
	    "or its sequence file in %s.\n",
	    target.name, target.time, target.pslTable, target.seqTable,
	    target.extFileTable, target.seqFile);
return NULL;
}

struct targetDb *targetDbLookup(char *db, char *name)
/* Given the name of a genomic database and the name of a PCR target
 * (or NULL to get all available PCR targets for db), query the
 * central database targetDb table and load the results.  Remove 
 * entries that are out of sync or have missing tables. */
{
struct sqlConnection *conn = hConnectCentral();
if (! sqlTableExists(conn, "targetDb"))
    {
    hDisconnectCentral(&conn);
    return NULL;
    }

struct targetDb *targetList = NULL;
struct sqlConnection *conn2 = hAllocConn(db);
struct sqlResult *sr;
char **row;
char query[2048];
safef(query, sizeof(query), "select * from targetDb where db = '%s' "
      "%s%s%s order by priority", db,
      isNotEmpty(name) ? "and name = '" : "",
      isNotEmpty(name) ? name : "",
      isNotEmpty(name) ? "' " : "");

sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    {
    struct targetDb *newTarg = targetDbMaybeLoad(conn2, row);
    if (newTarg)
	slAddHead(&targetList, newTarg);
    }
hFreeConn(&conn2);
hDisconnectCentral(&conn);
return targetList;
}

char *targetDbSetting(struct targetDb *tdb, char *name)
/* Return setting string or NULL if none exists. */
{
if (tdb == NULL)
    errAbort("Program error: null tdb passed to targetDbSetting.");
if (tdb->settingsHash == NULL)
    tdb->settingsHash = raFromString(tdb->settings);
return hashFindVal(tdb->settingsHash, name);
}

