/* parSpec.h was originally generated by the autoSql program, which also 
 * generated parSpec.c and parSpec.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef PARSPEC_H
#define PARSPEC_H

#define PARSPEC_NUM_COLS 7

struct parSpec
/* PAR specifications */
    {
    struct parSpec *next;  /* Next in singly linked list. */
    char *name;	/* PAR region name */
    char *chromA;	/* PAR region A chrom */
    int startA;	/* PAR region A start */
    int endA;	/* PAR region A end */
    char *chromB;	/* PAR region B chrom */
    int startB;	/* PAR region B start */
    int endB;	/* PAR region B end */
    };

struct parSpec *parSpecLoad(char **row);
/* Load a parSpec from row fetched with select * from parSpec
 * from database.  Dispose of this with parSpecFree(). */

struct parSpec *parSpecLoadAll(char *fileName);
/* Load all parSpec from whitespace-separated file.
 * Dispose of this with parSpecFreeList(). */

struct parSpec *parSpecLoadAllByChar(char *fileName, char chopper);
/* Load all parSpec from chopper separated file.
 * Dispose of this with parSpecFreeList(). */

#define parSpecLoadAllByTab(a) parSpecLoadAllByChar(a, '\t');
/* Load all parSpec from tab separated file.
 * Dispose of this with parSpecFreeList(). */

struct parSpec *parSpecCommaIn(char **pS, struct parSpec *ret);
/* Create a parSpec out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new parSpec */

void parSpecFree(struct parSpec **pEl);
/* Free a single dynamically allocated parSpec such as created
 * with parSpecLoad(). */

void parSpecFreeList(struct parSpec **pList);
/* Free a list of dynamically allocated parSpec's */

void parSpecOutput(struct parSpec *el, FILE *f, char sep, char lastSep);
/* Print out parSpec.  Separate fields with sep. Follow last field with lastSep. */

#define parSpecTabOut(el,f) parSpecOutput(el,f,'\t','\n');
/* Print out parSpec as a line in a tab-separated file. */

#define parSpecCommaOut(el,f) parSpecOutput(el,f,',',',');
/* Print out parSpec as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* PARSPEC_H */

