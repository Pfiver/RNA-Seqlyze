/* altGraph.h was originally generated by the autoSql program, which also 
 * generated altGraph.c and altGraph.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef ALTGRAPH_H
#define ALTGRAPH_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

struct altGraph
/* An alternatively spliced gene graph. */
    {
    struct altGraph *next;  /* Next in singly linked list. */
    unsigned id;	/* Unique ID */
    char *tName;	/* name of target sequence, often a chrom. */
    int tStart;	/* First bac touched by graph */
    int tEnd;	/* Start position in first bac */
    char strand[3];	/* + or - strand */
    unsigned vertexCount;	/* Number of vertices in graph */
    unsigned char *vTypes;	/* Type for each vertex */
    int *vPositions;	/* Position in target for each vertex */
    unsigned edgeCount;	/* Number of edges in graph */
    int *edgeStarts;	/* Array with start vertex of edges */
    int *edgeEnds;	/* Array with end vertex of edges. */
    int mrnaRefCount;	/* Number of supporting mRNAs. */
    char **mrnaRefs;	/* Ids of mrnas supporting this. */
    };

struct altGraph *altGraphLoad(char **row);
/* Load a altGraph from row fetched with select * from altGraph
 * from database.  Dispose of this with altGraphFree(). */

struct altGraph *altGraphLoadAll(char *fileName);
/* Load all altGraph from a tab-separated file.
 * Dispose of this with altGraphFreeList(). */

struct altGraph *altGraphLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all altGraph from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with altGraphFreeList(). */

void altGraphSaveToDb(struct sqlConnection *conn, struct altGraph *el, char *tableName, int updateSize);
/* Save altGraph as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use altGraphSaveToDbEscaped() */

void altGraphSaveToDbEscaped(struct sqlConnection *conn, struct altGraph *el, char *tableName, int updateSize);
/* Save altGraph as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than altGraphSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 

struct altGraph *altGraphCommaIn(char **pS, struct altGraph *ret);
/* Create a altGraph out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new altGraph */

void altGraphFree(struct altGraph **pEl);
/* Free a single dynamically allocated altGraph such as created
 * with altGraphLoad(). */

void altGraphFreeList(struct altGraph **pList);
/* Free a list of dynamically allocated altGraph's */

void altGraphOutput(struct altGraph *el, FILE *f, char sep, char lastSep);
/* Print out altGraph.  Separate fields with sep. Follow last field with lastSep. */

#define altGraphTabOut(el,f) altGraphOutput(el,f,'\t','\n');
/* Print out altGraph as a line in a tab-separated file. */

#define altGraphCommaOut(el,f) altGraphOutput(el,f,',',',');
/* Print out altGraph as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

int altGraphNumAltSplices(struct altGraph *ag);
/* Count number of times that exons have more than one edge through them */

#endif /* ALTGRAPH_H */

