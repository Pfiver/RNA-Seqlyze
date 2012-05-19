/* cnpLocke.h was originally generated by the autoSql program, which also 
 * generated cnpLocke.c and cnpLocke.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef CNPLOCKE_H
#define CNPLOCKE_H

#define CNPLOCKE_NUM_COLS 5

struct cnpLocke
/* CNP data from Locke */
    {
    struct cnpLocke *next;  /* Next in singly linked list. */
    char *chrom;	/* Reference sequence chromosome or scaffold */
    unsigned chromStart;	/* Start position in chrom */
    unsigned chromEnd;	/* End position in chrom */
    char *name;	/* BAC Name */
    char *variationType;	/* {Gain},{Loss},{Gain and Loss} */
    };

void cnpLockeStaticLoad(char **row, struct cnpLocke *ret);
/* Load a row from cnpLocke table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cnpLocke *cnpLockeLoad(char **row);
/* Load a cnpLocke from row fetched with select * from cnpLocke
 * from database.  Dispose of this with cnpLockeFree(). */

struct cnpLocke *cnpLockeLoadAll(char *fileName);
/* Load all cnpLocke from whitespace-separated file.
 * Dispose of this with cnpLockeFreeList(). */

struct cnpLocke *cnpLockeLoadAllByChar(char *fileName, char chopper);
/* Load all cnpLocke from chopper separated file.
 * Dispose of this with cnpLockeFreeList(). */

#define cnpLockeLoadAllByTab(a) cnpLockeLoadAllByChar(a, '\t');
/* Load all cnpLocke from tab separated file.
 * Dispose of this with cnpLockeFreeList(). */

struct cnpLocke *cnpLockeCommaIn(char **pS, struct cnpLocke *ret);
/* Create a cnpLocke out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cnpLocke */

void cnpLockeFree(struct cnpLocke **pEl);
/* Free a single dynamically allocated cnpLocke such as created
 * with cnpLockeLoad(). */

void cnpLockeFreeList(struct cnpLocke **pList);
/* Free a list of dynamically allocated cnpLocke's */

void cnpLockeOutput(struct cnpLocke *el, FILE *f, char sep, char lastSep);
/* Print out cnpLocke.  Separate fields with sep. Follow last field with lastSep. */

#define cnpLockeTabOut(el,f) cnpLockeOutput(el,f,'\t','\n');
/* Print out cnpLocke as a line in a tab-separated file. */

#define cnpLockeCommaOut(el,f) cnpLockeOutput(el,f,',',',');
/* Print out cnpLocke as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* CNPLOCKE_H */

