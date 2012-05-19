/* cnpSharpSample.c was originally generated by the autoSql program, which also 
 * generated cnpSharpSample.h and cnpSharpSample.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "cnpSharpSample.h"


void cnpSharpSampleStaticLoad(char **row, struct cnpSharpSample *ret)
/* Load a row from cnpSharpSample table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->bac = row[0];
ret->sample = row[1];
ret->batch = sqlUnsigned(row[2]);
ret->value = atof(row[3]);
ret->gender = row[4];
}

struct cnpSharpSample *cnpSharpSampleLoad(char **row)
/* Load a cnpSharpSample from row fetched with select * from cnpSharpSample
 * from database.  Dispose of this with cnpSharpSampleFree(). */
{
struct cnpSharpSample *ret;

AllocVar(ret);
ret->bac = cloneString(row[0]);
ret->sample = cloneString(row[1]);
ret->batch = sqlUnsigned(row[2]);
ret->value = atof(row[3]);
ret->gender = cloneString(row[4]);
return ret;
}

struct cnpSharpSample *cnpSharpSampleLoadAll(char *fileName) 
/* Load all cnpSharpSample from a whitespace-separated file.
 * Dispose of this with cnpSharpSampleFreeList(). */
{
struct cnpSharpSample *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[5];

while (lineFileRow(lf, row))
    {
    el = cnpSharpSampleLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct cnpSharpSample *cnpSharpSampleLoadAllByChar(char *fileName, char chopper) 
/* Load all cnpSharpSample from a chopper separated file.
 * Dispose of this with cnpSharpSampleFreeList(). */
{
struct cnpSharpSample *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[5];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = cnpSharpSampleLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct cnpSharpSample *cnpSharpSampleLoadByQuery(struct sqlConnection *conn, char *query)
/* Load all cnpSharpSample from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with cnpSharpSampleFreeList(). */
{
struct cnpSharpSample *list = NULL, *el;
struct sqlResult *sr;
char **row;

sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = cnpSharpSampleLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
return list;
}

void cnpSharpSampleSaveToDb(struct sqlConnection *conn, struct cnpSharpSample *el, char *tableName, int updateSize)
/* Save cnpSharpSample as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use cnpSharpSampleSaveToDbEscaped() */
{
struct dyString *update = newDyString(updateSize);
dyStringPrintf(update, "insert into %s values ( '%s','%s',%u,%g,'%s')", 
	tableName,  el->bac,  el->sample,  el->batch,  el->value,  el->gender);
sqlUpdate(conn, update->string);
freeDyString(&update);
}

void cnpSharpSampleSaveToDbEscaped(struct sqlConnection *conn, struct cnpSharpSample *el, char *tableName, int updateSize)
/* Save cnpSharpSample as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than cnpSharpSampleSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 
{
struct dyString *update = newDyString(updateSize);
char  *bac, *sample, *gender;
bac = sqlEscapeString(el->bac);
sample = sqlEscapeString(el->sample);
gender = sqlEscapeString(el->gender);

dyStringPrintf(update, "insert into %s values ( '%s','%s',%u,%g,'%s')", 
	tableName,  bac,  sample, el->batch , el->value ,  gender);
sqlUpdate(conn, update->string);
freeDyString(&update);
freez(&bac);
freez(&sample);
freez(&gender);
}

struct cnpSharpSample *cnpSharpSampleCommaIn(char **pS, struct cnpSharpSample *ret)
/* Create a cnpSharpSample out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cnpSharpSample */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->bac = sqlStringComma(&s);
ret->sample = sqlStringComma(&s);
ret->batch = sqlUnsignedComma(&s);
ret->value = sqlFloatComma(&s);
ret->gender = sqlStringComma(&s);
*pS = s;
return ret;
}

void cnpSharpSampleFree(struct cnpSharpSample **pEl)
/* Free a single dynamically allocated cnpSharpSample such as created
 * with cnpSharpSampleLoad(). */
{
struct cnpSharpSample *el;

if ((el = *pEl) == NULL) return;
freeMem(el->bac);
freeMem(el->sample);
freeMem(el->gender);
freez(pEl);
}

void cnpSharpSampleFreeList(struct cnpSharpSample **pList)
/* Free a list of dynamically allocated cnpSharpSample's */
{
struct cnpSharpSample *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    cnpSharpSampleFree(&el);
    }
*pList = NULL;
}

void cnpSharpSampleOutput(struct cnpSharpSample *el, FILE *f, char sep, char lastSep) 
/* Print out cnpSharpSample.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->bac);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->sample);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->batch);
fputc(sep,f);
fprintf(f, "%g", el->value);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->gender);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

