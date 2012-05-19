/* wikiTrack.c was originally generated by the autoSql program, which also 
 * generated wikiTrack.h and wikiTrack.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "wikiTrack.h"
#include "errCatch.h"


void wikiTrackStaticLoad(char **row, struct wikiTrack *ret)
/* Load a row from wikiTrack table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->bin = sqlUnsigned(row[0]);
ret->chrom = row[1];
ret->chromStart = sqlUnsigned(row[2]);
ret->chromEnd = sqlUnsigned(row[3]);
ret->name = row[4];
ret->score = sqlUnsigned(row[5]);
safecpy(ret->strand, sizeof(ret->strand), row[6]);
ret->db = row[7];
ret->owner = row[8];
ret->color = row[9];
ret->class = row[10];
ret->creationDate = row[11];
ret->lastModifiedDate = row[12];
ret->descriptionKey = row[13];
ret->id = sqlUnsigned(row[14]);
ret->geneSymbol = row[15];
}

struct wikiTrack *wikiTrackLoad(char **row)
/* Load a wikiTrack from row fetched with select * from wikiTrack
 * from database.  Dispose of this with wikiTrackFree(). */
{
struct wikiTrack *ret;

AllocVar(ret);
ret->bin = sqlUnsigned(row[0]);
ret->chrom = cloneString(row[1]);
ret->chromStart = sqlUnsigned(row[2]);
ret->chromEnd = sqlUnsigned(row[3]);
ret->name = cloneString(row[4]);
ret->score = sqlUnsigned(row[5]);
safecpy(ret->strand, sizeof(ret->strand), row[6]);
ret->db = cloneString(row[7]);
ret->owner = cloneString(row[8]);
ret->color = cloneString(row[9]);
ret->class = cloneString(row[10]);
ret->creationDate = cloneString(row[11]);
ret->lastModifiedDate = cloneString(row[12]);
ret->descriptionKey = cloneString(row[13]);
ret->id = sqlUnsigned(row[14]);
ret->geneSymbol = cloneString(row[15]);
return ret;
}

struct wikiTrack *wikiTrackLoadAll(char *fileName) 
/* Load all wikiTrack from a whitespace-separated file.
 * Dispose of this with wikiTrackFreeList(). */
{
struct wikiTrack *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[16];

while (lineFileRow(lf, row))
    {
    el = wikiTrackLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wikiTrack *wikiTrackLoadAllByChar(char *fileName, char chopper) 
/* Load all wikiTrack from a chopper separated file.
 * Dispose of this with wikiTrackFreeList(). */
{
struct wikiTrack *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[16];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = wikiTrackLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wikiTrack *wikiTrackLoadByQuery(struct sqlConnection *conn, char *query)
/* Load all wikiTrack from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with wikiTrackFreeList(). */
{
struct wikiTrack *list = NULL, *el;
struct sqlResult *sr;
char **row;

sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = wikiTrackLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
return list;
}

void wikiTrackSaveToDb(struct sqlConnection *conn, struct wikiTrack *el, char *tableName, int updateSize)
/* Save wikiTrack as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use wikiTrackSaveToDbEscaped() */
{
struct dyString *update = newDyString(updateSize);
dyStringPrintf(update, "insert into %s values ( %u,'%s',%u,%u,'%s',%u,'%s','%s','%s','%s','%s','%s','%s','%s',%u,'%s')", 
	tableName,  el->bin,  el->chrom,  el->chromStart,  el->chromEnd,  el->name,  el->score,  el->strand,  el->db,  el->owner,  el->color,  el->class,  el->creationDate,  el->lastModifiedDate,  el->descriptionKey,  el->id,  el->geneSymbol);
sqlUpdate(conn, update->string);
freeDyString(&update);
}

void wikiTrackSaveToDbEscaped(struct sqlConnection *conn, struct wikiTrack *el, char *tableName, int updateSize)
/* Save wikiTrack as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than wikiTrackSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 
{
struct dyString *update = newDyString(updateSize);
char  *chrom, *name, *strand, *db, *owner, *color, *class, *creationDate, *lastModifiedDate, *descriptionKey, *geneSymbol;
chrom = sqlEscapeString(el->chrom);
name = sqlEscapeString(el->name);
strand = sqlEscapeString(el->strand);
db = sqlEscapeString(el->db);
owner = sqlEscapeString(el->owner);
color = sqlEscapeString(el->color);
class = sqlEscapeString(el->class);
creationDate = sqlEscapeString(el->creationDate);
lastModifiedDate = sqlEscapeString(el->lastModifiedDate);
descriptionKey = sqlEscapeString(el->descriptionKey);
geneSymbol = sqlEscapeString(el->geneSymbol);

dyStringPrintf(update, "insert into %s values ( %u,'%s',%u,%u,'%s',%u,'%s','%s','%s','%s','%s','%s','%s','%s',%u,'%s')", 
	tableName, el->bin ,  chrom, el->chromStart , el->chromEnd ,  name, el->score ,  strand,  db,  owner,  color,  class,  creationDate,  lastModifiedDate,  descriptionKey, el->id ,  geneSymbol);
sqlUpdate(conn, update->string);
freeDyString(&update);
freez(&chrom);
freez(&name);
freez(&strand);
freez(&db);
freez(&owner);
freez(&color);
freez(&class);
freez(&creationDate);
freez(&lastModifiedDate);
freez(&descriptionKey);
freez(&geneSymbol);
}

struct wikiTrack *wikiTrackCommaIn(char **pS, struct wikiTrack *ret)
/* Create a wikiTrack out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wikiTrack */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->bin = sqlUnsignedComma(&s);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
ret->db = sqlStringComma(&s);
ret->owner = sqlStringComma(&s);
ret->color = sqlStringComma(&s);
ret->class = sqlStringComma(&s);
ret->creationDate = sqlStringComma(&s);
ret->lastModifiedDate = sqlStringComma(&s);
ret->descriptionKey = sqlStringComma(&s);
ret->id = sqlUnsignedComma(&s);
ret->geneSymbol = sqlStringComma(&s);
*pS = s;
return ret;
}

void wikiTrackFree(struct wikiTrack **pEl)
/* Free a single dynamically allocated wikiTrack such as created
 * with wikiTrackLoad(). */
{
struct wikiTrack *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->db);
freeMem(el->owner);
freeMem(el->color);
freeMem(el->class);
freeMem(el->creationDate);
freeMem(el->lastModifiedDate);
freeMem(el->descriptionKey);
freeMem(el->geneSymbol);
freez(pEl);
}

void wikiTrackFreeList(struct wikiTrack **pList)
/* Free a list of dynamically allocated wikiTrack's */
{
struct wikiTrack *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    wikiTrackFree(&el);
    }
*pList = NULL;
}

void wikiTrackOutput(struct wikiTrack *el, FILE *f, char sep, char lastSep) 
/* Print out wikiTrack.  Separate fields with sep. Follow last field with lastSep. */
{
fprintf(f, "%u", el->bin);
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
fprintf(f, "%s", el->db);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->owner);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->color);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->class);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->creationDate);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->lastModifiedDate);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->descriptionKey);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->id);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->geneSymbol);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#include "hgConfig.h"
#include "wikiLink.h"
#include "cheapcgi.h"
#include "hdb.h"
#include "hui.h"
#include "net.h"
#include "cart.h"
#include "hPrint.h"
#include "grp.h"
#include "obscure.h"
#include "hCommon.h"
#include "web.h"
#include "hgColors.h"

#ifdef NOT
static void savePosInTextBox(char *chrom, int start, int end)
/* Save basic position/database info in text box and hidden var. 
   Positions becomes chrom:start-end*/
{
char position[128];
char *newPos;
snprintf(position, 128, "%s:%d-%d", chrom, start, end);
newPos = addCommasToPos(position);
cgiMakeTextVar("getDnaPos", newPos, strlen(newPos) + 2);
cgiContinueHiddenVar("db");
}
#endif

boolean wikiTrackEnabled(char *database, char **wikiUserName)
/*determine if wikiTrack can be used, and is this user logged into the wiki ?*/
{
static boolean done = FALSE;
static boolean status = FALSE;
static boolean wikiUp = TRUE;
static char *userName = NULL;

/* do not repeat this query */
if (done)
    {
    if (wikiUserName)
	*wikiUserName = userName;  /* returning name indicates logged in */
    return status;
    }

done = TRUE;	/*	only need to do this once	*/

if (wikiUserName)
    *wikiUserName = NULL;  /* assume not logged in until proven otherwise */

/*	potentially a comma separated list */

char *dbListString = cfgOption(CFG_WIKI_DB_LIST);
char *dbList[256];
int dbCount = 0;
if (NULL != dbListString)
    {
    dbCount = chopByChar(cloneString(dbListString), ',', dbList,
	ArraySize(dbList));
    }
boolean validDb = TRUE;
if (dbCount > 0)
    {
    validDb = FALSE;
    int i;
    for (i = 0; i < dbCount; ++i)
	{
	if (sameWord(dbList[i],database))
	    {
	    validDb = TRUE;
	    break;
	    }
	}
    }
/* must have wiki login system enabled, and the new cfg options exist too. */
if (validDb && wikiLinkEnabled() &&
	(cfgOption(CFG_WIKI_SESSION_COOKIE) != NULL) &&
	(cfgOption(CFG_WIKI_BROWSER) != NULL) &&
	(cfgOption(CFG_WIKI_URL) != NULL))
    {
    char *wikiUser = wikiLinkUserName();
    if ( (wikiUser) &&
	(findCookieData(cfgOption(CFG_WIKI_SESSION_COOKIE)) != NULL) )
	{
	    struct htmlPage *page = NULL;
	    /* protect against temporarily offline wiki site */
	    struct errCatch *errCatch = errCatchNew();
	    if (errCatchStart(errCatch))
		{
		page = fetchEditPage(TEST_EMAIL_VERIFIED);
		}
	    errCatchEnd(errCatch);
	    if (errCatch->gotError) // we think it is supposed to be there
		{
		wikiUp = FALSE;		// but it will not respond
		}
	    errCatchFree(&errCatch);
	    char *loginExpired = NULL;
	    if (page)
		loginExpired = stringIn(LOGIN_EXPIRED, page->fullText);
	    if (loginExpired == NULL)
		userName = wikiUser;	/* save result for next time */
	}
    /* see if table exists, create it if it is not yet there */
    struct sqlConnection *wikiConn = wikiConnect();
    if (! sqlTableExists(wikiConn,WIKI_TRACK_TABLE))
	{
	char *query = wikiTrackGetCreateSql(WIKI_TRACK_TABLE);
	sqlUpdate(wikiConn, query);
	freeMem(query);
	}
    wikiDisconnect(&wikiConn);
    if (wikiUp)
	status = TRUE; /* system is enabled */
    }
if (wikiUserName)
    *wikiUserName = userName;  /* returning name indicates logged in */

return status;
}

/* from hg/lib/wikiTrack.sql */
static char *createString =
"CREATE TABLE %s (\n"
    "bin smallint unsigned not null,\n"
    "chrom varchar(255) not null,\n"
    "chromStart int unsigned not null,\n"
    "chromEnd int unsigned not null,\n"
    "name varchar(255) not null,\n"
    "score int unsigned not null,\n"
    "strand char(1) not null,\n"
    "db varchar(36) not null,\n"
    "owner varchar(255) not null,\n"
    "color varchar(255) not null,\n"
    "class varchar(255) not null,\n"
    "creationDate varchar(255) not null,\n"
    "lastModifiedDate varchar(255) not null,\n"
    "descriptionKey varchar(255) not null,\n"
    "id int unsigned not null auto_increment,\n"
    "geneSymbol varchar(255) null,\n"
    "PRIMARY KEY(id),\n"
    "INDEX chrom (db,bin,chrom),\n"
    "INDEX name (db,name),\n"
    "INDEX gene (geneSymbol)\n"
")\n";

char *wikiTrackGetCreateSql(char *tableName)
/* return sql create statement for wiki track with tableName */
{
struct dyString *createTable = dyStringNew(512);

dyStringPrintf(createTable, createString, tableName);

return (dyStringCannibalize(&createTable));
}

char *wikiDbName()
/* return name of database where wiki track is located
    currently this is central.db but the future may be configurable */
{
static char *dbName = NULL;

if (dbName)
    return dbName;

char setting[64];
safef(setting, sizeof(setting), "central.db");
dbName = cfgOption(setting);
return dbName;
}

struct sqlConnection *wikiConnect()
/* connect to db where wikiTrack table is located
 *	currently this is hConnectCentral() but the future may be
 *	configurable */
{
struct sqlConnection *conn = hConnectCentral();
return conn;
}

void wikiDisconnect(struct sqlConnection **pConn)
/* disconnect from wikiTrack table database */
{
hDisconnectCentral(pConn);
}

struct wikiTrack *findWikiItemId(char *wikiItemId)
/* given a wikiItemId return the row from the table */
{
struct wikiTrack *item;
char query[256];
struct sqlConnection *wikiConn = wikiConnect();

safef(query, ArraySize(query), "SELECT * FROM %s WHERE id='%s' limit 1",
	WIKI_TRACK_TABLE, wikiItemId);

item = wikiTrackLoadByQuery(wikiConn, query);
if (NULL == item)
    errAbort("display wiki item: failed to load item '%s'", wikiItemId);
wikiDisconnect(&wikiConn);

return item;
}

struct wikiTrack *findWikiItemByGeneSymbol(char *db, char *geneSymbol)
/* given a db and UCSC known gene geneSymbol, find the wiki item */
{
struct wikiTrack *item = NULL;

/* make sure neither of these arguments is NULL */
if (db && geneSymbol)
    {
    char query[256];
    struct sqlConnection *wikiConn = wikiConnect();
    safef(query, ArraySize(query),
	"SELECT * FROM %s WHERE db='%s' AND geneSymbol='%s' limit 1",
	    WIKI_TRACK_TABLE, db, geneSymbol);

    item = wikiTrackLoadByQuery(wikiConn, query);

    wikiDisconnect(&wikiConn);
    }

return item;
}

struct wikiTrack *findWikiItemByName(char *db, char *name)
/* given a db,name pair return the row from the table, can return NULL */
{
struct wikiTrack *item = NULL;

/* make sure neither of these arguments is NULL */
if (name && db)
    {
    char query[256];
    struct sqlConnection *wikiConn = wikiConnect();
    safef(query, ArraySize(query),
	"SELECT * FROM %s WHERE db='%s' AND name='%s' limit 1",
	    WIKI_TRACK_TABLE, db, name);

    item = wikiTrackLoadByQuery(wikiConn, query);

    wikiDisconnect(&wikiConn);
    }

return item;
}

static char *stripEditURLs(char *rendered)
/* test for actual text, remove edit sections and any html comment strings */
{
char *stripped = cloneString(rendered);
char *found = NULL;
char *begin = "<div class=\"editsection\"";
char *end = "></a>";

/* XXXX is this response going to be language dependent ? */
if (stringIn(WIKI_NO_TEXT_RESPONSE,rendered))
	return NULL;

/* remove any edit sections */
while (NULL != (found = stringBetween(begin, end, stripped)) )
    {
    if (strlen(found) > 0)
	{
	struct dyString *rm = newDyString(1024);
	dyStringPrintf(rm, "%s%s%s", begin, found, end);
	stripString(stripped, rm->string);
	freeMem(found);
	freeDyString(&rm);
	}
    }

/* and remove comment strings from the wiki */
begin = "<!--";
end = "-->";
while (NULL != (found = stringBetween(begin, end, stripped)) )
    {
    if (strlen(found) > 0)
	{
	struct dyString *rm = newDyString(1024);
	dyStringPrintf(rm, "%s%s%s", begin, found, end);
	stripString(stripped, rm->string);
	freeMem(found);
	freeDyString(&rm);
	}
    }

return stripped;
}

void htmlCloneFormVarSet(struct htmlForm *parent,
	struct htmlForm *clone, char *name, char *val)
/* clone form variable from parent, with new value,
 * if *val is NULL, clone val from parent
 */
{
struct htmlFormVar *cloneVar;
struct htmlFormVar *var;
if (parent == NULL)
    errAbort("Null parent form passed to htmlCloneFormVarSet");
if (clone == NULL)
    errAbort("Null clone form passed to htmlCloneFormVarSet");
var = htmlFormVarGet(parent, name);
if (var == NULL)
    errAbort("Variable '%s' not found in parent in htmlCloneFormVarSet", name);

AllocVar(cloneVar);
cloneVar->name = cloneString(var->name);
cloneVar->tagName = cloneString(var->tagName);
cloneVar->type = cloneString(var->type);
if (NULL == val)
    cloneVar->curVal = cloneString(var->curVal);
else
    cloneVar->curVal = cloneString(val);
slAddHead(&clone->vars, cloneVar);

}

char *fetchWikiRawText(char *descriptionKey)
/* fetch page from wiki in raw form as it is in the edit form */
{
char wikiPageUrl[512];
safef(wikiPageUrl, sizeof(wikiPageUrl), "%s/index.php?title=%s&action=raw",
	cfgOptionDefault(CFG_WIKI_URL, NULL), descriptionKey);
struct lineFile *lf = netLineFileMayOpen(wikiPageUrl);

struct dyString *wikiPage = newDyString(1024);
if (lf)
    {
    char *line;
    int lineSize;
    while (lineFileNext(lf, &line, &lineSize))
	dyStringPrintf(wikiPage, "%s\n", line);
    lineFileClose(&lf);
    }

/* test for text, remove any edit sections and comment strings */
char *rawText = NULL;
if (wikiPage->string)
    {
    /* XXXX is this response going to be language dependent ? */
    if (stringIn(WIKI_NO_TEXT_RESPONSE,wikiPage->string))
	return NULL;
    rawText = dyStringCannibalize(&wikiPage);
    }
freeDyString(&wikiPage);

return rawText;
}

char *fetchWikiRenderedText(char *descriptionKey)
/* fetch page from wiki in rendered form, strip it of edit URLs,
 *	html comments, and test for actual proper return.
 *  returned string can be freed after use */
{
char wikiPageUrl[512];
safef(wikiPageUrl, sizeof(wikiPageUrl), "%s/index.php/%s?action=render",
	cfgOptionDefault(CFG_WIKI_URL, NULL), descriptionKey);
struct lineFile *lf = netLineFileMayOpen(wikiPageUrl);

struct dyString *wikiPage = newDyString(1024);
if (lf)
    {
    char *line;
    int lineSize;
    while (lineFileNext(lf, &line, &lineSize))
	dyStringPrintf(wikiPage, "%s\n", line);
    lineFileClose(&lf);
    }

/* test for text, remove any edit sections and comment strings */
char *strippedRender = NULL;
if (wikiPage->string)
    strippedRender = stripEditURLs(wikiPage->string);
freeDyString(&wikiPage);

return strippedRender;
}

char *adjustWikiUrls (char *comments) 
/* change relative to full urls for images */
{
char *com;
char buf[128];
char *url = cfgOptionDefault(CFG_WIKI_URL, NULL);
safef(buf, sizeof(buf), "img src=\"%s/images", url);
com = replaceChars(comments, "img src=\"/images", buf);
return com;
}

void displayComments(struct wikiTrack *item)
/* display the rendered comments for this item */
{
char *url = cfgOptionDefault(CFG_WIKI_URL, NULL);

hPrintf("<B>Description and comments from the wiki article "
  "<A HREF=\"%s/index.php/%s\" TARGET=_blank>%s:</B></A><BR>\n",
       url, item->descriptionKey, item->descriptionKey);

/* the nbsp is displayed if the browser does not honor the frame tag */
hPrintf("<P>\n<IFRAME SRC='%s/index.php/%s?action=render'\n"
	" style='width:800px;height:600px;' frameborder='1'\n"
	" scrolling='yes' marginwidth='5' marginheight='5'\n"
	"&nbsp;</IFRAME>\n</P>\n", url, item->descriptionKey);

#ifdef NOT
char *comments = fetchWikiRenderedText(item->descriptionKey);
if (comments && (strlen(comments) > 2))
    {
    /* change relative to full urls for images */
    char *com = adjustWikiUrls(comments);
    hPrintf("\n%s\n", com);
    }
else
    hPrintf("\n(no comments for this item at the current time)<BR>\n");
#endif
}

struct htmlPage *fetchEditPage(char *descriptionKey)
/* fetch edit page for descriptionKey page name in wiki */
{
struct htmlCookie *cookie;
char wikiPageUrl[512];

/* must pass the session cookie from the wiki in order to edit */
AllocVar(cookie);
cookie->name = cloneString(cfgOption(CFG_WIKI_SESSION_COOKIE));
cookie->value = cloneString(findCookieData(cookie->name));

/* fetch the edit page to get the wpEditToken, and current contents */
safef(wikiPageUrl, sizeof(wikiPageUrl), "%s/index.php/%s?action=edit",
	cfgOptionDefault(CFG_WIKI_URL, NULL), descriptionKey);

char *fullText = htmlSlurpWithCookies(wikiPageUrl,cookie);
struct htmlPage *page = htmlPageParseOk(wikiPageUrl, fullText);
/* fullText pointer is placed in page->fullText */
/* the submit on the edit is going to need these cookies */
page->cookies = cookie;

return (page);
}

void prefixComments(struct wikiTrack *item, char *comments, char *userName,
    char *seqName, int winStart, int winEnd, char *database,
	char *extraHeader, char *extraTag, char *category)
/* add comments at the beginning of an existing wiki item */
{
/* do nothing if given nothing */
if (NULL == comments)
    return;

struct dyString *content = newDyString(1024);

struct htmlPage *page = fetchEditPage(item->descriptionKey);

/* create a stripped down edit form, we don't want all the variables */
struct htmlForm *strippedEditForm;
AllocVar(strippedEditForm);

struct htmlForm *currentEditForm = htmlFormGet(page,"editform");
if (NULL == currentEditForm)
    errAbort("prefixComments: can not get editform ? (wiki login confused ?)");

strippedEditForm->name = cloneString(currentEditForm->name);
/* the lower case "post" in the editform does not work ? */
/*strippedEditForm->method = cloneString(currentEditForm->method);*/
strippedEditForm->method = cloneString("POST");
strippedEditForm->startTag = currentEditForm->startTag;
strippedEditForm->endTag = currentEditForm->endTag;

/* fetch any current page contents in the edit form to continue them */
struct htmlFormVar *wpTextbox1 =
	htmlPageGetVar(page, currentEditForm, "wpTextbox1");

/* these new comments are the first thing on the page */
dyStringPrintf(content, "''comments added: ~~~~''");
if (extraTag)
    dyStringPrintf(content, "&nbsp;'''%s'''", extraTag);
dyStringPrintf(content, "\n\n%s<BR>\n", comments);

/* might want to recreate the header if it has gone missing */
/* decide on whether adding comments to existing text, or starting a
 *	new article from scratch.
 *	This function could be extended to actually checking the current
 *	contents to see if the "Category:" or "created:" lines have been
 *	removed, and then restore them.
 */
if (!(wpTextbox1->curVal && (strlen(wpTextbox1->curVal) > 2)))
    {
    boolean recreateHeader = FALSE;
    char position[128];
    char *newPos;
    char *userSignature;
    /* In the case where this is a restoration of the header lines,
     *	may be a different creator than this user adding comments.
     *	So, get the header line correct to represent the actual creator.
     */
    if (sameWord(userName, item->owner))
	userSignature = cloneString("~~~~");
    else
	{
	struct dyString *tt = newDyString(1024);
	dyStringPrintf(tt, "[[User:%s|%s]] ", item->owner, item->owner);
	dyStringPrintf(tt, "%s", item->creationDate);
	userSignature = dyStringCannibalize(&tt);
	recreateHeader = TRUE;
	}
    snprintf(position, 128, "%s:%d-%d", seqName, winStart+1, winEnd);
    newPos = addCommasToPos(database, position);
    
    if (extraHeader)
	{
	dyStringPrintf(content, "%s\n%s\n",
	    category, extraHeader);
	}
    else
	{
	dyStringPrintf(content, "%s\n"
"[http://%s/cgi-bin/hgTracks?db=%s&wikiTrack=pack&position=%s:%d-%d %s %s]"
	"&nbsp;&nbsp;<B>'%s'</B>&nbsp;&nbsp;",
	category,
	    cfgOptionDefault(CFG_WIKI_BROWSER, DEFAULT_BROWSER), database,
		seqName, winStart+1, winEnd, database, newPos, item->name);
	}
    dyStringPrintf(content, "''created: %s''", userSignature);
    if (extraTag)
	dyStringPrintf(content, "&nbsp;'''%s'''", extraTag);
    }

/* append previous content, if any */
if (wpTextbox1->curVal && (strlen(wpTextbox1->curVal) > 2))
    {
    char *rawText = fetchWikiRawText(item->descriptionKey);
    dyStringPrintf(content, "\n\n%s", rawText);
    }
else
    {
    dyStringPrintf(content, "\n\n%s", NO_ITEM_COMMENT_SUPPLIED);
    }

htmlCloneFormVarSet(currentEditForm, strippedEditForm,
	"wpTextbox1", content->string);
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSummary", "");
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSection", "");
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpMinoredit", "1");
/*
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSave", "Save page");
*/
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpEdittime", NULL);
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpEditToken", NULL);

htmlPageSetVar(page,currentEditForm, "wpTextbox1", content->string);
htmlPageSetVar(page,currentEditForm, "wpSummary", "");
htmlPageSetVar(page,currentEditForm, "wpSection", "");
htmlPageSetVar(page,currentEditForm, "wpMinoredit", "1");
htmlPageSetVar(page,currentEditForm, "wpSave", "Save page");

char newUrl[1024];
char *wikiHost = wikiLinkHost();
/* fake out htmlPageFromForm since it doesn't understand the colon : */
safef(newUrl, ArraySize(newUrl), "http://%s%s",
	wikiHost, currentEditForm->action);
freeMem(wikiHost);
/* something, somewhere encoded the & into &amp; which does not work */
char *fixedString = replaceChars(newUrl, "&amp;", "&");
currentEditForm->action = cloneString(fixedString);
strippedEditForm->action = cloneString(fixedString);
struct htmlPage *editPage = htmlPageFromForm(page,strippedEditForm,"submit", "Submit");
freeMem(fixedString);

if (NULL == editPage)
    errAbort("prefixComments: the edit is failing ?");

freeDyString(&content);
}	/*	void prefixComments()	*/


void addDescription(struct wikiTrack *item, char *userName,
    char *seqName, int winStart, int winEnd, struct cart *cart,
	char *database, char *extraHeader, char *extraTag, char *category)
/* add description to the end of an existing wiki item */
{
char *newComments = cartNonemptyString(cart, NEW_ITEM_COMMENT);
struct dyString *content = newDyString(1024);

/* was nothing changed in the add comments entry box ? */
if (sameWord(ADD_ITEM_COMMENT_DEFAULT,newComments))
    return;

struct htmlPage *page = fetchEditPage(item->descriptionKey);

/* create a stripped down edit form, we don't want all the variables */
struct htmlForm *strippedEditForm;
AllocVar(strippedEditForm);

struct htmlForm *currentEditForm = htmlFormGet(page,"editform");
if (NULL == currentEditForm)
    errAbort("addDescription: can not get editform ? (wiki login confused ?)");

strippedEditForm->name = cloneString(currentEditForm->name);
/* the lower case "post" in the editform does not work ? */
/*strippedEditForm->method = cloneString(currentEditForm->method);*/
strippedEditForm->method = cloneString("POST");
strippedEditForm->startTag = currentEditForm->startTag;
strippedEditForm->endTag = currentEditForm->endTag;

/* fetch any current page contents in the edit form to continue them */
struct htmlFormVar *wpTextbox1 =
	htmlPageGetVar(page, currentEditForm, "wpTextbox1");

/* decide on whether adding comments to existing text, or starting a
 *	new article from scratch.
 *	This function could be extended to actually checking the current
 *	contents to see if the "Category:" or "created:" lines have been
 *	removed, and then restore them.
 */
if (wpTextbox1->curVal && (strlen(wpTextbox1->curVal) > 2))
    {
    char *rawText = fetchWikiRawText(item->descriptionKey);
    dyStringPrintf(content, "%s\n\n''comments added: ~~~~''", rawText);
    if (extraTag)
	dyStringPrintf(content, "&nbsp;'''%s'''\n\n", extraTag);
    else
	dyStringPrintf(content, "\n\n");
    }
else
    {
    boolean recreateHeader = FALSE;
    char position[128];
    char *newPos;
    char *userSignature;
    /* In the case where this is a restoration of the header lines,
     *	may be a different creator than this user adding comments.
     *	So, get the header line correct to represent the actual creator.
     */
    if (sameWord(userName, item->owner))
	userSignature = cloneString("~~~~");
    else
	{
	struct dyString *tt = newDyString(1024);
	dyStringPrintf(tt, "[[User:%s|%s]] ", item->owner, item->owner);
	dyStringPrintf(tt, "%s", item->creationDate);
	userSignature = dyStringCannibalize(&tt);
	recreateHeader = TRUE;
	}
    snprintf(position, 128, "%s:%d-%d", seqName, winStart+1, winEnd);
    newPos = addCommasToPos(database, position);
    if (extraHeader)
	{
	dyStringPrintf(content, "%s\n%s\n",
	    category, extraHeader);
	}
    else
	{
	dyStringPrintf(content, "%s\n"
"[http://%s/cgi-bin/hgTracks?db=%s&wikiTrack=pack&position=%s:%d-%d %s %s]"
	"&nbsp;&nbsp;<B>'%s'</B>&nbsp;&nbsp;",
	category,
	    cfgOptionDefault(CFG_WIKI_BROWSER, DEFAULT_BROWSER), database,
		seqName, winStart+1, winEnd, database, newPos, item->name);
	}
    dyStringPrintf(content, "''created: %s''", userSignature);
    if (extraTag)
	dyStringPrintf(content, "&nbsp;'''%s'''", extraTag);
    dyStringPrintf(content, "\n\n");
    if (recreateHeader)
	{
	dyStringPrintf(content, "''comments added: ~~~~''");
	if (extraTag)
	    dyStringPrintf(content, "&nbsp;'''%s'''\n\n", extraTag);
	else
	    dyStringPrintf(content, "\n\n");
	}
    }

if (sameWord(NEW_ITEM_COMMENT_DEFAULT,newComments))
    dyStringPrintf(content, "%s\n\n", NO_ITEM_COMMENT_SUPPLIED);
else
    dyStringPrintf(content, "%s\n\n", newComments);

htmlCloneFormVarSet(currentEditForm, strippedEditForm,
	"wpTextbox1", content->string);
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSummary", "");
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSection", "");
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpMinoredit", "1");
/*
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpSave", "Save page");
*/
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpEdittime", NULL);
htmlCloneFormVarSet(currentEditForm, strippedEditForm, "wpEditToken", NULL);

htmlPageSetVar(page,currentEditForm, "wpTextbox1", content->string);
htmlPageSetVar(page,currentEditForm, "wpSummary", "");
htmlPageSetVar(page,currentEditForm, "wpSection", "");
htmlPageSetVar(page,currentEditForm, "wpMinoredit", "1");
htmlPageSetVar(page,currentEditForm, "wpSave", "Save page");

char newUrl[1024];
char *wikiHost = wikiLinkHost();
/* fake out htmlPageFromForm since it doesn't understand the colon : */
safef(newUrl, ArraySize(newUrl), "http://%s%s",
	wikiHost, currentEditForm->action);
freeMem(wikiHost);
/* something, somewhere encoded the & into &amp; which does not work */
char *fixedString = replaceChars(newUrl, "&amp;", "&");
currentEditForm->action = cloneString(fixedString);
strippedEditForm->action = cloneString(fixedString);
struct htmlPage *editPage = htmlPageFromForm(page,strippedEditForm,"submit", "Submit");
freeMem(fixedString);

if (NULL == editPage)
    errAbort("addDescription: the edit is failing ?");

freeDyString(&content);
}	/*	void addDescription()	*/

char *encodedReturnUrl(char *(*hgUrl)())
/* Return a CGI-encoded URL with hgsid.  Free when done.
 *	The given function hgUrl() will construct the actual cgi binary URL
 */
{
char retBuf[1024];
safef(retBuf, sizeof(retBuf), "http://%s/%s", cgiServerName(), (*hgUrl)());
return cgiEncode(retBuf);
}   

boolean emailVerified(boolean showMessage)
/* TRUE indicates email has been verified for this wiki user */
{
struct htmlPage *page = fetchEditPage(TEST_EMAIL_VERIFIED);
char *verifyEmail = stringIn(EMAIL_NEEDS_TO_BE_VERIFIED, page->fullText);
/* sometimes the genome browser thinks the user is logged in, but the
 *	wiki doesn't think so.  So, verify that an editform exists
 *	within the edit page to confirm edit is OK.
 */
struct htmlForm *currentEditForm = htmlFormGet(page,"editform");
htmlPageFree(&page);
if ((currentEditForm != NULL) && (verifyEmail == NULL))
    return TRUE;
else
    {
    if (showMessage)
	{
	hPrintf("<P>%s annotations.  %s "
	    "<A HREF=\"%s/index.php/Special:Preferences\" "
	    "TARGET=_blank>user preferences.</A></P>\n",
		EMAIL_NEEDS_TO_BE_VERIFIED, USER_PREFERENCES_MESSAGE, 
		    cfgOptionDefault(CFG_WIKI_URL, NULL));
	}
    return FALSE;
    }
}

boolean isWikiEditor(char *userName)
/* check if user name is on list of editors */
{
boolean editor = FALSE;

if (userName)
    {
    char *editors = cfgOptionDefault(CFG_WIKI_EDITORS, NULL);
    if (editors)
	{
	char *cloneEditors = cloneString(editors);
	int i;
	int wordCount = chopByChar(cloneEditors, ',', NULL, 0);
	char **words = (char **)needMem((size_t)(wordCount * sizeof(char *)));
	chopByChar(cloneEditors, ',', words, wordCount);
	for (i = 0; i < wordCount; ++i)
	    {
	    if (sameWord(userName, words[i]))
		{
		editor = TRUE;
		break;
		}
	    }
	freeMem(cloneEditors);
	}
    }
return editor;
} 

char *wikiUrl(struct wikiTrack *item)
/* construct a URL to the wiki page for the specified item
	free the returned string when done with it.  */
{
char *siteUrl = cfgOptionDefault(CFG_WIKI_URL, NULL);
struct dyString *itemUrl = dyStringNew(64);
dyStringPrintf(itemUrl, "%s/index.php/%s TARGET=_blank", siteUrl,
        item->descriptionKey);
return dyStringCannibalize(&itemUrl);
}
