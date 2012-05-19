/* taxonXref.c was originally generated by the autoSql program, which also 
 * generated taxonXref.h and taxonXref.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "taxonXref.h"


void taxonXrefStaticLoad(char **row, struct taxonXref *ret)
/* Load a row from taxonXref table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->organism = row[0];
ret->taxon = sqlUnsigned(row[1]);
ret->name = row[2];
ret->toGenus = row[3];
}

struct taxonXref *taxonXrefLoad(char **row)
/* Load a taxonXref from row fetched with select * from taxonXref
 * from database.  Dispose of this with taxonXrefFree(). */
{
struct taxonXref *ret;

AllocVar(ret);
ret->organism = cloneString(row[0]);
ret->taxon = sqlUnsigned(row[1]);
ret->name = cloneString(row[2]);
ret->toGenus = cloneString(row[3]);
return ret;
}

struct taxonXref *taxonXrefLoadAll(char *fileName) 
/* Load all taxonXref from a whitespace-separated file.
 * Dispose of this with taxonXrefFreeList(). */
{
struct taxonXref *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[4];

while (lineFileRow(lf, row))
    {
    el = taxonXrefLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct taxonXref *taxonXrefLoadAllByChar(char *fileName, char chopper) 
/* Load all taxonXref from a chopper separated file.
 * Dispose of this with taxonXrefFreeList(). */
{
struct taxonXref *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[4];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = taxonXrefLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct taxonXref *taxonXrefLoadByQuery(struct sqlConnection *conn, char *query)
/* Load all taxonXref from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with taxonXrefFreeList(). */
{
struct taxonXref *list = NULL, *el;
struct sqlResult *sr;
char **row;

sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = taxonXrefLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
return list;
}

void taxonXrefSaveToDb(struct sqlConnection *conn, struct taxonXref *el, char *tableName, int updateSize)
/* Save taxonXref as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use taxonXrefSaveToDbEscaped() */
{
struct dyString *update = newDyString(updateSize);
dyStringPrintf(update, "insert into %s values ( '%s',%u,'%s','%s')", 
	tableName,  el->organism,  el->taxon,  el->name,  el->toGenus);
sqlUpdate(conn, update->string);
freeDyString(&update);
}

void taxonXrefSaveToDbEscaped(struct sqlConnection *conn, struct taxonXref *el, char *tableName, int updateSize)
/* Save taxonXref as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than taxonXrefSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 
{
struct dyString *update = newDyString(updateSize);
char  *organism, *name, *toGenus;
organism = sqlEscapeString(el->organism);
name = sqlEscapeString(el->name);
toGenus = sqlEscapeString(el->toGenus);

dyStringPrintf(update, "insert into %s values ( '%s',%u,'%s','%s')", 
	tableName,  organism, el->taxon ,  name,  toGenus);
sqlUpdate(conn, update->string);
freeDyString(&update);
freez(&organism);
freez(&name);
freez(&toGenus);
}

struct taxonXref *taxonXrefCommaIn(char **pS, struct taxonXref *ret)
/* Create a taxonXref out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new taxonXref */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->organism = sqlStringComma(&s);
ret->taxon = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->toGenus = sqlStringComma(&s);
*pS = s;
return ret;
}

void taxonXrefFree(struct taxonXref **pEl)
/* Free a single dynamically allocated taxonXref such as created
 * with taxonXrefLoad(). */
{
struct taxonXref *el;

if ((el = *pEl) == NULL) return;
freeMem(el->organism);
freeMem(el->name);
freeMem(el->toGenus);
freez(pEl);
}

void taxonXrefFreeList(struct taxonXref **pList)
/* Free a list of dynamically allocated taxonXref's */
{
struct taxonXref *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    taxonXrefFree(&el);
    }
*pList = NULL;
}

void taxonXrefOutput(struct taxonXref *el, FILE *f, char sep, char lastSep) 
/* Print out taxonXref.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->organism);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->taxon);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->toGenus);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

