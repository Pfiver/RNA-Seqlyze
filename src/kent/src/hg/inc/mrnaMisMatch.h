/* mrnaMisMatch.h was originally generated by the autoSql program, which also 
 * generated mrnaMisMatch.c and mrnaMisMatch.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef MRNAMISMATCH_H
#define MRNAMISMATCH_H

#define MRNAMISMATCH_NUM_COLS 11

struct mrnaMisMatch
/* list of mismatches between mrna and multiple hits in the genome */
    {
    struct mrnaMisMatch *next;  /* Next in singly linked list. */
    char *name;	/* mRNA accession */
    char mrnaBase;	/* base in mrna */
    int mrnaLoc;	/* location of mismatch in mrna */
    int misMatchCount;	/* count of misMatches  */
    char *bases;	/* genomic bases for each misMatch */
    char **chroms;	/* chrom of each misMatch */
    unsigned *tStarts;	/* chrom start of each misMatch */
    char *strands;	/* strand of each misMatch */
    unsigned *loci;	/* loci index for each mismatch */
    int snpCount;	/* number of snps in this position */
    char **snps;	/* list of dbSNP ids */
    };

struct mrnaMisMatch *mrnaMisMatchLoad(char **row);
/* Load a mrnaMisMatch from row fetched with select * from mrnaMisMatch
 * from database.  Dispose of this with mrnaMisMatchFree(). */

struct mrnaMisMatch *mrnaMisMatchLoadAll(char *fileName);
/* Load all mrnaMisMatch from whitespace-separated file.
 * Dispose of this with mrnaMisMatchFreeList(). */

struct mrnaMisMatch *mrnaMisMatchLoadAllByChar(char *fileName, char chopper);
/* Load all mrnaMisMatch from chopper separated file.
 * Dispose of this with mrnaMisMatchFreeList(). */

#define mrnaMisMatchLoadAllByTab(a) mrnaMisMatchLoadAllByChar(a, '\t');
/* Load all mrnaMisMatch from tab separated file.
 * Dispose of this with mrnaMisMatchFreeList(). */

struct mrnaMisMatch *mrnaMisMatchCommaIn(char **pS, struct mrnaMisMatch *ret);
/* Create a mrnaMisMatch out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new mrnaMisMatch */

void mrnaMisMatchFree(struct mrnaMisMatch **pEl);
/* Free a single dynamically allocated mrnaMisMatch such as created
 * with mrnaMisMatchLoad(). */

void mrnaMisMatchFreeList(struct mrnaMisMatch **pList);
/* Free a list of dynamically allocated mrnaMisMatch's */

void mrnaMisMatchOutput(struct mrnaMisMatch *el, FILE *f, char sep, char lastSep);
/* Print out mrnaMisMatch.  Separate fields with sep. Follow last field with lastSep. */

#define mrnaMisMatchTabOut(el,f) mrnaMisMatchOutput(el,f,'\t','\n');
/* Print out mrnaMisMatch as a line in a tab-separated file. */

#define mrnaMisMatchCommaOut(el,f) mrnaMisMatchOutput(el,f,',',',');
/* Print out mrnaMisMatch as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* MRNAMISMATCH_H */

