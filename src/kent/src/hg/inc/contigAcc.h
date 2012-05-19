/* contigAcc.h was originally generated by the autoSql program, which also 
 * generated contigAcc.c and contigAcc.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef CONTIGACC_H
#define CONTIGACC_H

#define CONTIGACC_NUM_COLS 2

struct contigAcc
/* Map a contig to its accession. */
    {
    struct contigAcc *next;  /* Next in singly linked list. */
    char *contig;	/* Contig name */
    char *acc;	/* Genbank accession */
    };

void contigAccStaticLoad(char **row, struct contigAcc *ret);
/* Load a row from contigAcc table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct contigAcc *contigAccLoad(char **row);
/* Load a contigAcc from row fetched with select * from contigAcc
 * from database.  Dispose of this with contigAccFree(). */

struct contigAcc *contigAccLoadAll(char *fileName);
/* Load all contigAcc from whitespace-separated file.
 * Dispose of this with contigAccFreeList(). */

struct contigAcc *contigAccLoadAllByChar(char *fileName, char chopper);
/* Load all contigAcc from chopper separated file.
 * Dispose of this with contigAccFreeList(). */

#define contigAccLoadAllByTab(a) contigAccLoadAllByChar(a, '\t');
/* Load all contigAcc from tab separated file.
 * Dispose of this with contigAccFreeList(). */

struct contigAcc *contigAccCommaIn(char **pS, struct contigAcc *ret);
/* Create a contigAcc out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new contigAcc */

void contigAccFree(struct contigAcc **pEl);
/* Free a single dynamically allocated contigAcc such as created
 * with contigAccLoad(). */

void contigAccFreeList(struct contigAcc **pList);
/* Free a list of dynamically allocated contigAcc's */

void contigAccOutput(struct contigAcc *el, FILE *f, char sep, char lastSep);
/* Print out contigAcc.  Separate fields with sep. Follow last field with lastSep. */

#define contigAccTabOut(el,f) contigAccOutput(el,f,'\t','\n');
/* Print out contigAcc as a line in a tab-separated file. */

#define contigAccCommaOut(el,f) contigAccOutput(el,f,',',',');
/* Print out contigAcc as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* CONTIGACC_H */

