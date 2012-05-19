/* wabAli.h was originally generated by the autoSql program, which also 
 * generated wabAli.c and wabAli.sql.  This header links the database and the RAM 
 * representation of objects. */

#ifndef WABALI_H
#define WABALI_H

struct wabAli
/* Information on a waba alignment */
    {
    struct wabAli *next;  /* Next in singly linked list. */
    char *query;	/* Foreign sequence name. */
    unsigned qStart;	/* Start in query (other species) */
    unsigned qEnd;	/* End in query. */
    char strand[2];	/* + or - relative orientation. */
    char *chrom;	/* Chromosome (current species). */
    unsigned chromStart;	/* Start in chromosome. */
    unsigned chromEnd;	/* End in chromosome. */
    unsigned milliScore;	/* Base identity in parts per thousand. */
    unsigned symCount;	/* Number of symbols in alignment. */
    char *qSym;	/* Query bases plus '-'s. */
    char *tSym;	/* Target bases plus '-'s. */
    char *hSym;	/* Hidden Markov symbols. */
    };

void wabAliStaticLoad(char **row, struct wabAli *ret);
/* Load a row from wabAli table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct wabAli *wabAliLoad(char **row);
/* Load a wabAli from row fetched with select * from wabAli
 * from database.  Dispose of this with wabAliFree(). */

struct wabAli *wabAliCommaIn(char **pS, struct wabAli *ret);
/* Create a wabAli out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wabAli */

void wabAliFree(struct wabAli **pEl);
/* Free a single dynamically allocated wabAli such as created
 * with wabAliLoad(). */

void wabAliFreeList(struct wabAli **pList);
/* Free a list of dynamically allocated wabAli's */

void wabAliOutput(struct wabAli *el, FILE *f, char sep, char lastSep);
/* Print out wabAli.  Separate fields with sep. Follow last field with lastSep. */

#define wabAliTabOut(el,f) wabAliOutput(el,f,'\t','\n');
/* Print out wabAli as a line in a tab-separated file. */

#define wabAliCommaOut(el,f) wabAliOutput(el,f,',',',');
/* Print out wabAli as a comma separated list including final comma. */

#endif /* WABALI_H */

