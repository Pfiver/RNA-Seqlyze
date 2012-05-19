/* tigrCmrGene.c was originally generated by the autoSql program, which also 
 * generated tigrCmrGene.h and tigrCmrGene.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "tigrCmrGene.h"


void tigrCmrGeneStaticLoad(char **row, struct tigrCmrGene *ret)
/* Load a row from tigrCmrGene table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->bin = sqlSigned(row[0]);
ret->chrom = row[1];
ret->chromStart = sqlUnsigned(row[2]);
ret->chromEnd = sqlUnsigned(row[3]);
ret->name = row[4];
ret->score = sqlUnsigned(row[5]);
strcpy(ret->strand, row[6]);
ret->tigrCommon = row[7];
ret->tigrGene = row[8];
ret->tigrECN = row[9];
ret->primLocus = row[10];
ret->tigrLength = sqlUnsigned(row[11]);
ret->tigrPepLength = sqlUnsigned(row[12]);
ret->tigrMainRole = row[13];
ret->tigrSubRole = row[14];
ret->swissProt = row[15];
ret->genbank = row[16];
ret->tigrMw = atof(row[17]);
ret->tigrPi = atof(row[18]);
ret->tigrGc = atof(row[19]);
ret->goTerm = row[20];
}

struct tigrCmrGene *tigrCmrGeneLoad(char **row)
/* Load a tigrCmrGene from row fetched with select * from tigrCmrGene
 * from database.  Dispose of this with tigrCmrGeneFree(). */
{
struct tigrCmrGene *ret;

AllocVar(ret);
ret->bin = sqlSigned(row[0]);
ret->chrom = cloneString(row[1]);
ret->chromStart = sqlUnsigned(row[2]);
ret->chromEnd = sqlUnsigned(row[3]);
ret->name = cloneString(row[4]);
ret->score = sqlUnsigned(row[5]);
strcpy(ret->strand, row[6]);
ret->tigrCommon = cloneString(row[7]);
ret->tigrGene = cloneString(row[8]);
ret->tigrECN = cloneString(row[9]);
ret->primLocus = cloneString(row[10]);
ret->tigrLength = sqlUnsigned(row[11]);
ret->tigrPepLength = sqlUnsigned(row[12]);
ret->tigrMainRole = cloneString(row[13]);
ret->tigrSubRole = cloneString(row[14]);
ret->swissProt = cloneString(row[15]);
ret->genbank = cloneString(row[16]);
ret->tigrMw = atof(row[17]);
ret->tigrPi = atof(row[18]);
ret->tigrGc = atof(row[19]);
ret->goTerm = cloneString(row[20]);
return ret;
}

struct tigrCmrGene *tigrCmrGeneLoadAll(char *fileName) 
/* Load all tigrCmrGene from a whitespace-separated file.
 * Dispose of this with tigrCmrGeneFreeList(). */
{
struct tigrCmrGene *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[21];

while (lineFileRow(lf, row))
    {
    el = tigrCmrGeneLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct tigrCmrGene *tigrCmrGeneLoadAllByChar(char *fileName, char chopper) 
/* Load all tigrCmrGene from a chopper separated file.
 * Dispose of this with tigrCmrGeneFreeList(). */
{
struct tigrCmrGene *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[21];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = tigrCmrGeneLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct tigrCmrGene *tigrCmrGeneLoadByQuery(struct sqlConnection *conn, char *query)
/* Load all tigrCmrGene from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with tigrCmrGeneFreeList(). */
{
struct tigrCmrGene *list = NULL, *el;
struct sqlResult *sr;
char **row;

sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = tigrCmrGeneLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
return list;
}

void tigrCmrGeneSaveToDb(struct sqlConnection *conn, struct tigrCmrGene *el, char *tableName, int updateSize)
/* Save tigrCmrGene as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use tigrCmrGeneSaveToDbEscaped() */
{
struct dyString *update = newDyString(updateSize);
dyStringPrintf(update, "insert into %s values ( %d,'%s',%u,%u,'%s',%u,'%s',%s,'%s','%s','%s',%u,%u,%s,%s,'%s','%s',%f,%f,%f,'%s')", 
	tableName,  el->bin,  el->chrom,  el->chromStart,  el->chromEnd,  el->name,  el->score,  el->strand,  el->tigrCommon,  el->tigrGene,  el->tigrECN,  el->primLocus,  el->tigrLength,  el->tigrPepLength,  el->tigrMainRole,  el->tigrSubRole,  el->swissProt,  el->genbank,  el->tigrMw,  el->tigrPi,  el->tigrGc,  el->goTerm);
sqlUpdate(conn, update->string);
freeDyString(&update);
}

void tigrCmrGeneSaveToDbEscaped(struct sqlConnection *conn, struct tigrCmrGene *el, char *tableName, int updateSize)
/* Save tigrCmrGene as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than tigrCmrGeneSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 
{
struct dyString *update = newDyString(updateSize);
char  *chrom, *name, *strand, *tigrCommon, *tigrGene, *tigrECN, *primLocus, *tigrMainRole, *tigrSubRole, *swissProt, *genbank, *goTerm;
chrom = sqlEscapeString(el->chrom);
name = sqlEscapeString(el->name);
strand = sqlEscapeString(el->strand);
tigrCommon = sqlEscapeString(el->tigrCommon);
tigrGene = sqlEscapeString(el->tigrGene);
tigrECN = sqlEscapeString(el->tigrECN);
primLocus = sqlEscapeString(el->primLocus);
tigrMainRole = sqlEscapeString(el->tigrMainRole);
tigrSubRole = sqlEscapeString(el->tigrSubRole);
swissProt = sqlEscapeString(el->swissProt);
genbank = sqlEscapeString(el->genbank);
goTerm = sqlEscapeString(el->goTerm);

dyStringPrintf(update, "insert into %s values ( %d,'%s',%u,%u,'%s',%u,'%s','%s','%s','%s','%s',%u,%u,'%s','%s','%s','%s',%f,%f,%f,'%s')", 
	tableName, el->bin ,  chrom, el->chromStart , el->chromEnd ,  name, el->score ,  strand,  tigrCommon,  tigrGene,  tigrECN,  primLocus, el->tigrLength , el->tigrPepLength ,  tigrMainRole,  tigrSubRole,  swissProt,  genbank, el->tigrMw , el->tigrPi , el->tigrGc ,  goTerm);
sqlUpdate(conn, update->string);
freeDyString(&update);
freez(&chrom);
freez(&name);
freez(&strand);
freez(&tigrCommon);
freez(&tigrGene);
freez(&tigrECN);
freez(&primLocus);
freez(&tigrMainRole);
freez(&tigrSubRole);
freez(&swissProt);
freez(&genbank);
freez(&goTerm);
}

struct tigrCmrGene *tigrCmrGeneCommaIn(char **pS, struct tigrCmrGene *ret)
/* Create a tigrCmrGene out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new tigrCmrGene */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->bin = sqlSignedComma(&s);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
ret->tigrCommon = sqlStringComma(&s);
ret->tigrGene = sqlStringComma(&s);
ret->tigrECN = sqlStringComma(&s);
ret->primLocus = sqlStringComma(&s);
ret->tigrLength = sqlUnsignedComma(&s);
ret->tigrPepLength = sqlUnsignedComma(&s);
ret->tigrMainRole = sqlStringComma(&s);
ret->tigrSubRole = sqlStringComma(&s);
ret->swissProt = sqlStringComma(&s);
ret->genbank = sqlStringComma(&s);
ret->tigrMw = sqlFloatComma(&s);
ret->tigrPi = sqlFloatComma(&s);
ret->tigrGc = sqlFloatComma(&s);
ret->goTerm = sqlStringComma(&s);
*pS = s;
return ret;
}

void tigrCmrGeneFree(struct tigrCmrGene **pEl)
/* Free a single dynamically allocated tigrCmrGene such as created
 * with tigrCmrGeneLoad(). */
{
struct tigrCmrGene *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->tigrCommon);
freeMem(el->tigrGene);
freeMem(el->tigrECN);
freeMem(el->primLocus);
freeMem(el->tigrMainRole);
freeMem(el->tigrSubRole);
freeMem(el->swissProt);
freeMem(el->genbank);
freeMem(el->goTerm);
freez(pEl);
}

void tigrCmrGeneFreeList(struct tigrCmrGene **pList)
/* Free a list of dynamically allocated tigrCmrGene's */
{
struct tigrCmrGene *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    tigrCmrGeneFree(&el);
    }
*pList = NULL;
}

void tigrCmrGeneOutput(struct tigrCmrGene *el, FILE *f, char sep, char lastSep) 
/* Print out tigrCmrGene.  Separate fields with sep. Follow last field with lastSep. */
{
fprintf(f, "%d", el->bin);
fputc(sep,f);
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
fprintf(f, "%u", el->score);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->tigrCommon);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->tigrGene);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->tigrECN);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->primLocus);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->tigrLength);
fputc(sep,f);
fprintf(f, "%u", el->tigrPepLength);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->tigrMainRole);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->tigrSubRole);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->swissProt);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->genbank);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%f", el->tigrMw);
fputc(sep,f);
fprintf(f, "%f", el->tigrPi);
fputc(sep,f);
fprintf(f, "%f", el->tigrGc);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->goTerm);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

