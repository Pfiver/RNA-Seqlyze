/* cdsOrtho.h was originally generated by the autoSql program, which also 
 * generated cdsOrtho.c and cdsOrtho.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef CDSORTHO_H
#define CDSORTHO_H

#define CDSORTHO_NUM_COLS 8

struct cdsOrtho
/* Information about a CDS region in another species, created by looking at multiple alignment. */
    {
    struct cdsOrtho *next;  /* Next in singly linked list. */
    char *name;	/* Name of transcript */
    int start;	/* CDS start within transcript */
    int end;	/* CDS end within transcript */
    char *species;	/* Other species (or species database) */
    int missing;	/* Number of bases missing (non-aligning) */
    int orthoSize;	/* Size of orf in other species */
    int possibleSize;	/* Possible size of orf in other species */
    double ratio;	/* orthoSize/possibleSize */
    };

void cdsOrthoStaticLoad(char **row, struct cdsOrtho *ret);
/* Load a row from cdsOrtho table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cdsOrtho *cdsOrthoLoad(char **row);
/* Load a cdsOrtho from row fetched with select * from cdsOrtho
 * from database.  Dispose of this with cdsOrthoFree(). */

struct cdsOrtho *cdsOrthoLoadAll(char *fileName);
/* Load all cdsOrtho from whitespace-separated file.
 * Dispose of this with cdsOrthoFreeList(). */

struct cdsOrtho *cdsOrthoLoadAllByChar(char *fileName, char chopper);
/* Load all cdsOrtho from chopper separated file.
 * Dispose of this with cdsOrthoFreeList(). */

#define cdsOrthoLoadAllByTab(a) cdsOrthoLoadAllByChar(a, '\t');
/* Load all cdsOrtho from tab separated file.
 * Dispose of this with cdsOrthoFreeList(). */

struct cdsOrtho *cdsOrthoCommaIn(char **pS, struct cdsOrtho *ret);
/* Create a cdsOrtho out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdsOrtho */

void cdsOrthoFree(struct cdsOrtho **pEl);
/* Free a single dynamically allocated cdsOrtho such as created
 * with cdsOrthoLoad(). */

void cdsOrthoFreeList(struct cdsOrtho **pList);
/* Free a list of dynamically allocated cdsOrtho's */

void cdsOrthoOutput(struct cdsOrtho *el, FILE *f, char sep, char lastSep);
/* Print out cdsOrtho.  Separate fields with sep. Follow last field with lastSep. */

#define cdsOrthoTabOut(el,f) cdsOrthoOutput(el,f,'\t','\n');
/* Print out cdsOrtho as a line in a tab-separated file. */

#define cdsOrthoCommaOut(el,f) cdsOrthoOutput(el,f,',',',');
/* Print out cdsOrtho as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* CDSORTHO_H */

