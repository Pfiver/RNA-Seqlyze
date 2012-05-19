/* encodeErgeStableTransf.h was originally generated by the autoSql program, which also 
 * generated encodeErgeStableTransf.c and encodeErgeStableTransf.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef ENCODEERGESTABLETRANSF_H
#define ENCODEERGESTABLETRANSF_H

#define ENCODEERGESTABLETRANSF_NUM_COLS 14

struct encodeErgeStableTransf
/* ENCODE experimental data from dbERGEII */
    {
    struct encodeErgeStableTransf *next;  /* Next in singly linked list. */
    char *chrom;	/* Human chromosome */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    char *name;	/* Name of read - up to 255 characters */
    unsigned score;	/* Score from 0-1000.  1000 is best */
    char strand[2];	/* Value should be + or - */
    unsigned thickStart;	/* Start of where display should be thick (start codon) */
    unsigned thickEnd;	/* End of where display should be thick (stop codon) */
    unsigned reserved;	/* Always zero for now */
    unsigned blockCount;	/* Number of separate blocks (regions without gaps) */
    unsigned *blockSizes;	/* Comma separated list of block sizes */
    unsigned *chromStarts;	/* Start position of each block in relative to chromStart */
    char *Id;	/* dbERGEII Id */
    char *color;	/* RGB color values */
    };

struct encodeErgeStableTransf *encodeErgeStableTransfLoad(char **row);
/* Load a encodeErgeStableTransf from row fetched with select * from encodeErgeStableTransf
 * from database.  Dispose of this with encodeErgeStableTransfFree(). */

struct encodeErgeStableTransf *encodeErgeStableTransfLoadAll(char *fileName);
/* Load all encodeErgeStableTransf from whitespace-separated file.
 * Dispose of this with encodeErgeStableTransfFreeList(). */

struct encodeErgeStableTransf *encodeErgeStableTransfLoadAllByChar(char *fileName, char chopper);
/* Load all encodeErgeStableTransf from chopper separated file.
 * Dispose of this with encodeErgeStableTransfFreeList(). */

#define encodeErgeStableTransfLoadAllByTab(a) encodeErgeStableTransfLoadAllByChar(a, '\t');
/* Load all encodeErgeStableTransf from tab separated file.
 * Dispose of this with encodeErgeStableTransfFreeList(). */

struct encodeErgeStableTransf *encodeErgeStableTransfCommaIn(char **pS, struct encodeErgeStableTransf *ret);
/* Create a encodeErgeStableTransf out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new encodeErgeStableTransf */

void encodeErgeStableTransfFree(struct encodeErgeStableTransf **pEl);
/* Free a single dynamically allocated encodeErgeStableTransf such as created
 * with encodeErgeStableTransfLoad(). */

void encodeErgeStableTransfFreeList(struct encodeErgeStableTransf **pList);
/* Free a list of dynamically allocated encodeErgeStableTransf's */

void encodeErgeStableTransfOutput(struct encodeErgeStableTransf *el, FILE *f, char sep, char lastSep);
/* Print out encodeErgeStableTransf.  Separate fields with sep. Follow last field with lastSep. */

#define encodeErgeStableTransfTabOut(el,f) encodeErgeStableTransfOutput(el,f,'\t','\n');
/* Print out encodeErgeStableTransf as a line in a tab-separated file. */

#define encodeErgeStableTransfCommaOut(el,f) encodeErgeStableTransfOutput(el,f,',',',');
/* Print out encodeErgeStableTransf as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* ENCODEERGESTABLETRANSF_H */

