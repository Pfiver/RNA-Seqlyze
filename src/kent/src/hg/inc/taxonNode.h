/* taxonNode.h was originally generated by the autoSql program, which also 
 * generated taxonNode.c and taxonNode.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef TAXONNODE_H
#define TAXONNODE_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

#define TAXONNODE_NUM_COLS 13

struct taxonNode
/* ncbi taxonomy node tree */
    {
    struct taxonNode *next;  /* Next in singly linked list. */
    unsigned taxon;	/* node id in GenBank taxonomy database */
    unsigned parent;	/* parent node id in GenBank taxonomy database */
    char *rank;	/* rank of this node (superkingdom, kingdom, ...)  */
    char *emblcode;	/* locus-name prefix; not unique */
    unsigned division;	/* ncbiDivision id (0=Bacteria, 2=Mammal, 5=Primate, 6=Rodent, 10=Vertabrate...) */
    unsigned inheritedDiv;	/* 1 if node inherits division from parent */
    unsigned geneticCode;	/* genetic code used by species, see ncbiGencode */
    unsigned inheritedGC;	/* 1 if node inherits genetic code from parent */
    unsigned mitoGeneticCode;	/* genetic code of mitochondria see ncbiGencode */
    unsigned inheritedMitoGC;	/* 1 if node inherits mitochondrial gencode from parent */
    unsigned GenBankHidden;	/* 1 if name is suppressed in GenBank entry lineage */
    unsigned notSequenced;	/* 1 if this subtree has no sequence data yet */
    char *comments;	/* free-text comments and citations */
    };

void taxonNodeStaticLoad(char **row, struct taxonNode *ret);
/* Load a row from taxonNode table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct taxonNode *taxonNodeLoad(char **row);
/* Load a taxonNode from row fetched with select * from taxonNode
 * from database.  Dispose of this with taxonNodeFree(). */

struct taxonNode *taxonNodeLoadAll(char *fileName);
/* Load all taxonNode from whitespace-separated file.
 * Dispose of this with taxonNodeFreeList(). */

struct taxonNode *taxonNodeLoadAllByChar(char *fileName, char chopper);
/* Load all taxonNode from chopper separated file.
 * Dispose of this with taxonNodeFreeList(). */

#define taxonNodeLoadAllByTab(a) taxonNodeLoadAllByChar(a, '\t');
/* Load all taxonNode from tab separated file.
 * Dispose of this with taxonNodeFreeList(). */

struct taxonNode *taxonNodeLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all taxonNode from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with taxonNodeFreeList(). */

void taxonNodeSaveToDb(struct sqlConnection *conn, struct taxonNode *el, char *tableName, int updateSize);
/* Save taxonNode as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include" 
 * If worried about this use taxonNodeSaveToDbEscaped() */

void taxonNodeSaveToDbEscaped(struct sqlConnection *conn, struct taxonNode *el, char *tableName, int updateSize);
/* Save taxonNode as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically 
 * escapes all simple strings (not arrays of string) but may be slower than taxonNodeSaveToDb().
 * For example automatically copies and converts: 
 * "autosql's features include" --> "autosql\'s features include" 
 * before inserting into database. */ 

struct taxonNode *taxonNodeCommaIn(char **pS, struct taxonNode *ret);
/* Create a taxonNode out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new taxonNode */

void taxonNodeFree(struct taxonNode **pEl);
/* Free a single dynamically allocated taxonNode such as created
 * with taxonNodeLoad(). */

void taxonNodeFreeList(struct taxonNode **pList);
/* Free a list of dynamically allocated taxonNode's */

void taxonNodeOutput(struct taxonNode *el, FILE *f, char sep, char lastSep);
/* Print out taxonNode.  Separate fields with sep. Follow last field with lastSep. */

#define taxonNodeTabOut(el,f) taxonNodeOutput(el,f,'\t','\n');
/* Print out taxonNode as a line in a tab-separated file. */

#define taxonNodeCommaOut(el,f) taxonNodeOutput(el,f,',',',');
/* Print out taxonNode as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* TAXONNODE_H */

