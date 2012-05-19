/* expData.c was originally generated by the autoSql program, which also 
 * generated expData.h and expData.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "expData.h"


struct expData *expDataLoad(char **row)
/* Load a expData from row fetched with select * from expData
 * from database.  Dispose of this with expDataFree(). */
{
struct expData *ret;
int sizeOne;

AllocVar(ret);
ret->expCount = sqlUnsigned(row[1]);
ret->name = cloneString(row[0]);
sqlFloatDynamicArray(row[2], &ret->expScores, &sizeOne);
assert(sizeOne == ret->expCount);
return ret;
}

struct expData *expDataLoadAll(char *fileName) 
/* Load all expData from a whitespace-separated file.
 * Dispose of this with expDataFreeList(). */
{
struct expData *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileRow(lf, row))
    {
    el = expDataLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct expData *expDataLoadAllByChar(char *fileName, char chopper) 
/* Load all expData from a chopper separated file.
 * Dispose of this with expDataFreeList(). */
{
struct expData *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = expDataLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct expData *expDataCommaIn(char **pS, struct expData *ret)
/* Create a expData out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new expData */
{
char *s = *pS;
int i;

if (ret == NULL)
    AllocVar(ret);
ret->name = sqlStringComma(&s);
ret->expCount = sqlUnsignedComma(&s);
s = sqlEatChar(s, '{');
AllocArray(ret->expScores, ret->expCount);
for (i=0; i<ret->expCount; ++i)
    {
    ret->expScores[i] = sqlFloatComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
*pS = s;
return ret;
}

void expDataFree(struct expData **pEl)
/* Free a single dynamically allocated expData such as created
 * with expDataLoad(). */
{
struct expData *el;

if ((el = *pEl) == NULL) return;
freeMem(el->name);
freeMem(el->expScores);
freez(pEl);
}

void expDataFreeList(struct expData **pList)
/* Free a list of dynamically allocated expData's */
{
struct expData *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    expDataFree(&el);
    }
*pList = NULL;
}

void expDataOutput(struct expData *el, FILE *f, char sep, char lastSep) 
/* Print out expData.  Separate fields with sep. Follow last field with lastSep. */
{
int i;
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->expCount);
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->expCount; ++i)
    {
    fprintf(f, "%0.3f", el->expScores[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

void expDataCreateTable(struct sqlConnection *conn, char *table)
/* Create table with given name. */
{
char query[512];

safef(query, sizeof(query),
"CREATE TABLE %s (\n"
"    name varchar(255) not null,\n"
"    expCount int unsigned not null,\n"
"    expScores longblob not null,\n"
"    INDEX(name(10))\n"
")\n",   table);
sqlRemakeTable(conn, table, query);
}

struct expData *expDataLoadTableLimit(struct sqlConnection *conn, char *table, int limitRows)
/* Same as expDataLoadTable, but limit to only loading limitRows # of rows. */
{
char query[256];
char **row;
int numLoaded = 0;
struct expData *exps = NULL;
struct sqlResult *sr = NULL;
if (limitRows < 0)
    return NULL;
safef(query, sizeof(query), "select name, expCount, expScores from %s", table);
sr = sqlGetResult(conn, query);
if (limitRows > 0)
    {
    while (((row = sqlNextRow(sr)) != NULL) && (numLoaded < limitRows))
	{
	struct expData *addMe = expDataLoad(row);
	slAddHead(&exps, addMe);
	numLoaded++;
	}
    }
else
    {
    while ((row = sqlNextRow(sr)) != NULL)
	{
	struct expData *addMe = expDataLoad(row);
	slAddHead(&exps, addMe);
	}
    }
slReverse(&exps);
sqlFreeResult(&sr);
return exps;
}

struct expData *expDataLoadTable(struct sqlConnection *conn, char *table)
/* Load all the rows of an SQL table (already connected to the database) */
/* into a list and return it. This should work on BED 15 tables as well */
/* as native expData tables. */
{
return expDataLoadTableLimit(conn, table, 0);
}

struct expData *expDataConnectAndLoadTable(char *database, char *table)
/* Same thing as expDataLoadTableConn, but it does the extra step of */
/* connecting to a database first. */
{
struct sqlConnection *conn = sqlConnect(database);
struct expData *exps = expDataLoadTable(conn, table);
sqlDisconnect(&conn);
return exps;
}
