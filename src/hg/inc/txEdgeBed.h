/* txEdgeBed.h was originally generated by the autoSql program, which also 
 * generated txEdgeBed.c and txEdgeBed.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef TXEDGEBED_H
#define TXEDGEBED_H

#define TXEDGEBED_NUM_COLS 9

enum txEdgeBedType
    {
    txEdgeBedExon = 0,
    txEdgeBedIntron = 1,
    };
struct txEdgeBed
/* A transcription edge with first fields bed format. */
    {
    struct txEdgeBed *next;  /* Next in singly linked list. */
    char *chrom;	/* Chromosome or contig name */
    int chromStart;	/* Start position, zero-based */
    int chromEnd;	/* End position, non-inclusive */
    char *name;	/* Name of evidence supporting edge */
    int score;	/* Score - 0-1000 */
    char strand[2];	/* Strand - either plus or minus */
    char startType[2];	/* [ or ( for hard or soft */
    enum txEdgeBedType type;	/* edge type */
    char endType[2];	/* ] or ) for hard or soft */
    };

void txEdgeBedStaticLoad(char **row, struct txEdgeBed *ret);
/* Load a row from txEdgeBed table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct txEdgeBed *txEdgeBedLoad(char **row);
/* Load a txEdgeBed from row fetched with select * from txEdgeBed
 * from database.  Dispose of this with txEdgeBedFree(). */

struct txEdgeBed *txEdgeBedLoadAll(char *fileName);
/* Load all txEdgeBed from whitespace-separated file.
 * Dispose of this with txEdgeBedFreeList(). */

struct txEdgeBed *txEdgeBedLoadAllByChar(char *fileName, char chopper);
/* Load all txEdgeBed from chopper separated file.
 * Dispose of this with txEdgeBedFreeList(). */

#define txEdgeBedLoadAllByTab(a) txEdgeBedLoadAllByChar(a, '\t');
/* Load all txEdgeBed from tab separated file.
 * Dispose of this with txEdgeBedFreeList(). */

struct txEdgeBed *txEdgeBedCommaIn(char **pS, struct txEdgeBed *ret);
/* Create a txEdgeBed out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new txEdgeBed */

void txEdgeBedFree(struct txEdgeBed **pEl);
/* Free a single dynamically allocated txEdgeBed such as created
 * with txEdgeBedLoad(). */

void txEdgeBedFreeList(struct txEdgeBed **pList);
/* Free a list of dynamically allocated txEdgeBed's */

void txEdgeBedOutput(struct txEdgeBed *el, FILE *f, char sep, char lastSep);
/* Print out txEdgeBed.  Separate fields with sep. Follow last field with lastSep. */

#define txEdgeBedTabOut(el,f) txEdgeBedOutput(el,f,'\t','\n');
/* Print out txEdgeBed as a line in a tab-separated file. */

#define txEdgeBedCommaOut(el,f) txEdgeBedOutput(el,f,',',',');
/* Print out txEdgeBed as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* TXEDGEBED_H */

