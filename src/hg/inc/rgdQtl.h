/* rgdQtl.h was originally generated by the autoSql program, which also 
 * generated rgdQtl.c and rgdQtl.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef RGDQTL_H
#define RGDQTL_H

#define RGDQTL_NUM_COLS 5

struct rgdQtl
/* RGD QTL */
    {
    struct rgdQtl *next;  /* Next in singly linked list. */
    short bin;	/* bin for browser speed up */
    char *chrom;	/* Reference sequence chromosome or scaffold */
    unsigned chromStart;	/* Start in chromosome */
    unsigned chromEnd;	/* End in chromosome */
    char *name;	/* Name of QTL */
    };

void rgdQtlStaticLoad(char **row, struct rgdQtl *ret);
/* Load a row from rgdQtl table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct rgdQtl *rgdQtlLoad(char **row);
/* Load a rgdQtl from row fetched with select * from rgdQtl
 * from database.  Dispose of this with rgdQtlFree(). */

struct rgdQtl *rgdQtlLoadAll(char *fileName);
/* Load all rgdQtl from whitespace-separated file.
 * Dispose of this with rgdQtlFreeList(). */

struct rgdQtl *rgdQtlLoadAllByChar(char *fileName, char chopper);
/* Load all rgdQtl from chopper separated file.
 * Dispose of this with rgdQtlFreeList(). */

#define rgdQtlLoadAllByTab(a) rgdQtlLoadAllByChar(a, '\t');
/* Load all rgdQtl from tab separated file.
 * Dispose of this with rgdQtlFreeList(). */

struct rgdQtl *rgdQtlCommaIn(char **pS, struct rgdQtl *ret);
/* Create a rgdQtl out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new rgdQtl */

void rgdQtlFree(struct rgdQtl **pEl);
/* Free a single dynamically allocated rgdQtl such as created
 * with rgdQtlLoad(). */

void rgdQtlFreeList(struct rgdQtl **pList);
/* Free a list of dynamically allocated rgdQtl's */

void rgdQtlOutput(struct rgdQtl *el, FILE *f, char sep, char lastSep);
/* Print out rgdQtl.  Separate fields with sep. Follow last field with lastSep. */

#define rgdQtlTabOut(el,f) rgdQtlOutput(el,f,'\t','\n');
/* Print out rgdQtl as a line in a tab-separated file. */

#define rgdQtlCommaOut(el,f) rgdQtlOutput(el,f,',',',');
/* Print out rgdQtl as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* RGDQTL_H */

