/* sangerGeneToWBGeneID.h was originally generated by the autoSql program, which also 
 * generated sangerGeneToWBGeneID.c and sangerGeneToWBGeneID.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef SANGERGENETOWBGENEID_H
#define SANGERGENETOWBGENEID_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

#define SANGERGENETOWBGENEID_NUM_COLS 2

struct sangerGeneToWBGeneID
/* sanger gene names to WormBase Gene IDs translation */
    {
    struct sangerGeneToWBGeneID *next;  /* Next in singly linked list. */
    char *sangerGene;	/* sangerGene name */
    char *WBGeneID;	/* WormBase Gene ID */
    };

void sangerGeneToWBGeneIDStaticLoad(char **row, struct sangerGeneToWBGeneID *ret);
/* Load a row from sangerGeneToWBGeneID table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct sangerGeneToWBGeneID *sangerGeneToWBGeneIDLoad(char **row);
/* Load a sangerGeneToWBGeneID from row fetched with select * from sangerGeneToWBGeneID
 * from database.  Dispose of this with sangerGeneToWBGeneIDFree(). */

struct sangerGeneToWBGeneID *sangerGeneToWBGeneIDLoadAll(char *fileName);
/* Load all sangerGeneToWBGeneID from whitespace-separated file.
 * Dispose of this with sangerGeneToWBGeneIDFreeList(). */

struct sangerGeneToWBGeneID *sangerGeneToWBGeneIDLoadAllByChar(char *fileName, char chopper);
/* Load all sangerGeneToWBGeneID from chopper separated file.
 * Dispose of this with sangerGeneToWBGeneIDFreeList(). */

#define sangerGeneToWBGeneIDLoadAllByTab(a) sangerGeneToWBGeneIDLoadAllByChar(a, '\t');
/* Load all sangerGeneToWBGeneID from tab separated file.
 * Dispose of this with sangerGeneToWBGeneIDFreeList(). */

struct sangerGeneToWBGeneID *sangerGeneToWBGeneIDLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all sangerGeneToWBGeneID from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with sangerGeneToWBGeneIDFreeList(). */

void sangerGeneToWBGeneIDSaveToDb(struct sqlConnection *conn, struct sangerGeneToWBGeneID *el, char *tableName, int updateSize);
/* Save sangerGeneToWBGeneID as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use sangerGeneToWBGeneIDSaveToDbEscaped() */

void sangerGeneToWBGeneIDSaveToDbEscaped(struct sqlConnection *conn, struct sangerGeneToWBGeneID *el, char *tableName, int updateSize);
/* Save sangerGeneToWBGeneID as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than sangerGeneToWBGeneIDSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 

struct sangerGeneToWBGeneID *sangerGeneToWBGeneIDCommaIn(char **pS, struct sangerGeneToWBGeneID *ret);
/* Create a sangerGeneToWBGeneID out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new sangerGeneToWBGeneID */

void sangerGeneToWBGeneIDFree(struct sangerGeneToWBGeneID **pEl);
/* Free a single dynamically allocated sangerGeneToWBGeneID such as created
 * with sangerGeneToWBGeneIDLoad(). */

void sangerGeneToWBGeneIDFreeList(struct sangerGeneToWBGeneID **pList);
/* Free a list of dynamically allocated sangerGeneToWBGeneID's */

void sangerGeneToWBGeneIDOutput(struct sangerGeneToWBGeneID *el, FILE *f, char sep, char lastSep);
/* Print out sangerGeneToWBGeneID.  Separate fields with sep. Follow last field with lastSep. */

#define sangerGeneToWBGeneIDTabOut(el,f) sangerGeneToWBGeneIDOutput(el,f,'\t','\n');
/* Print out sangerGeneToWBGeneID as a line in a tab-separated file. */

#define sangerGeneToWBGeneIDCommaOut(el,f) sangerGeneToWBGeneIDOutput(el,f,',',',');
/* Print out sangerGeneToWBGeneID as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* SANGERGENETOWBGENEID_H */

