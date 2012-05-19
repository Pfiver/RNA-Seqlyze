/* easyGene.h was originally generated by the autoSql program, which also 
 * generated easyGene.c and easyGene.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef EASYGENE_H
#define EASYGENE_H

#define EASYGENE_NUM_COLS 15

struct easyGene
/* BED-style EasyGene format */
    {
    struct easyGene *next;  /* Next in singly linked list. */
    char *chrom;	/* Chromosome or FPC contig */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    char *name;	/* Name of item */
    unsigned score;	/* Score */
    char strand[2];	/* Strand */
    char *feat;	/* Feature identifier */
    float R;	/* R value */
    unsigned frame;	/* Frame */
    char *orf;	/* ORF identifier */
    char *startCodon;	/* Start Codon */
    float logOdds;	/* Log odds */
    char *descriptor;	/* EasyGene descriptor */
    char *swissProt;	/* Swiss-Prot match */
    char genbank[2];	/* Y or N this is the genbank one */
    };

void easyGeneStaticLoad(char **row, struct easyGene *ret);
/* Load a row from easyGene table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct easyGene *easyGeneLoad(char **row);
/* Load a easyGene from row fetched with select * from easyGene
 * from database.  Dispose of this with easyGeneFree(). */

struct easyGene *easyGeneLoadAll(char *fileName);
/* Load all easyGene from whitespace-separated file.
 * Dispose of this with easyGeneFreeList(). */

struct easyGene *easyGeneLoadAllByChar(char *fileName, char chopper);
/* Load all easyGene from chopper separated file.
 * Dispose of this with easyGeneFreeList(). */

#define easyGeneLoadAllByTab(a) easyGeneLoadAllByChar(a, '\t');
/* Load all easyGene from tab separated file.
 * Dispose of this with easyGeneFreeList(). */

struct easyGene *easyGeneCommaIn(char **pS, struct easyGene *ret);
/* Create a easyGene out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new easyGene */

void easyGeneFree(struct easyGene **pEl);
/* Free a single dynamically allocated easyGene such as created
 * with easyGeneLoad(). */

void easyGeneFreeList(struct easyGene **pList);
/* Free a list of dynamically allocated easyGene's */

void easyGeneOutput(struct easyGene *el, FILE *f, char sep, char lastSep);
/* Print out easyGene.  Separate fields with sep. Follow last field with lastSep. */

#define easyGeneTabOut(el,f) easyGeneOutput(el,f,'\t','\n');
/* Print out easyGene as a line in a tab-separated file. */

#define easyGeneCommaOut(el,f) easyGeneOutput(el,f,',',',');
/* Print out easyGene as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* EASYGENE_H */

