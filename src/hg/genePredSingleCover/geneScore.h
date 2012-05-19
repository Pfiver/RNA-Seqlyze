/* geneScore.h was originally generated by the autoSql program, which also 
 * generated geneScore.c and geneScore.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef GENESCORE_H
#define GENESCORE_H

#define GENESCORE_NUM_COLS 4

struct geneScore
/* object for loading gene score files */
    {
    struct geneScore *next;  /* Next in singly linked list. */
    char *name;	/* gene name */
    char *chrom;	/* chromosome name */
    int txStart;	/* gene txStart */
    float score;	/* score for gene */
    };

void geneScoreStaticLoad(char **row, struct geneScore *ret);
/* Load a row from geneScore table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct geneScore *geneScoreLoad(char **row);
/* Load a geneScore from row fetched with select * from geneScore
 * from database.  Dispose of this with geneScoreFree(). */

struct geneScore *geneScoreLoadAll(char *fileName);
/* Load all geneScore from whitespace-separated file.
 * Dispose of this with geneScoreFreeList(). */

struct geneScore *geneScoreLoadAllByChar(char *fileName, char chopper);
/* Load all geneScore from chopper separated file.
 * Dispose of this with geneScoreFreeList(). */

#define geneScoreLoadAllByTab(a) geneScoreLoadAllByChar(a, '\t');
/* Load all geneScore from tab separated file.
 * Dispose of this with geneScoreFreeList(). */

struct geneScore *geneScoreCommaIn(char **pS, struct geneScore *ret);
/* Create a geneScore out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new geneScore */

void geneScoreFree(struct geneScore **pEl);
/* Free a single dynamically allocated geneScore such as created
 * with geneScoreLoad(). */

void geneScoreFreeList(struct geneScore **pList);
/* Free a list of dynamically allocated geneScore's */

void geneScoreOutput(struct geneScore *el, FILE *f, char sep, char lastSep);
/* Print out geneScore.  Separate fields with sep. Follow last field with lastSep. */

#define geneScoreTabOut(el,f) geneScoreOutput(el,f,'\t','\n');
/* Print out geneScore as a line in a tab-separated file. */

#define geneScoreCommaOut(el,f) geneScoreOutput(el,f,',',',');
/* Print out geneScore as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* GENESCORE_H */

