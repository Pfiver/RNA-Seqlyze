/* transMapGene.h was originally generated by the autoSql program, which also 
 * generated transMapGene.c and transMapGene.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef TRANSMAPGENE_H
#define TRANSMAPGENE_H

#define TRANSMAPGENE_NUM_COLS 4

struct transMapGene
/* shared, gene-specific transMap information.  This is also a cdsSpec object */
    {
    struct transMapGene *next;  /* Next in singly linked list. */
    char *id;	/* unique sequence id */
    char *cds;	/* CDS specification, in NCBI format. */
    char db[17];	/* source db */
    char *geneName;	/* gene name */
    };

void transMapGeneStaticLoad(char **row, struct transMapGene *ret);
/* Load a row from transMapGene table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct transMapGene *transMapGeneLoad(char **row);
/* Load a transMapGene from row fetched with select * from transMapGene
 * from database.  Dispose of this with transMapGeneFree(). */

struct transMapGene *transMapGeneLoadAll(char *fileName);
/* Load all transMapGene from whitespace-separated file.
 * Dispose of this with transMapGeneFreeList(). */

struct transMapGene *transMapGeneLoadAllByChar(char *fileName, char chopper);
/* Load all transMapGene from chopper separated file.
 * Dispose of this with transMapGeneFreeList(). */

#define transMapGeneLoadAllByTab(a) transMapGeneLoadAllByChar(a, '\t');
/* Load all transMapGene from tab separated file.
 * Dispose of this with transMapGeneFreeList(). */

struct transMapGene *transMapGeneCommaIn(char **pS, struct transMapGene *ret);
/* Create a transMapGene out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new transMapGene */

void transMapGeneFree(struct transMapGene **pEl);
/* Free a single dynamically allocated transMapGene such as created
 * with transMapGeneLoad(). */

void transMapGeneFreeList(struct transMapGene **pList);
/* Free a list of dynamically allocated transMapGene's */

void transMapGeneOutput(struct transMapGene *el, FILE *f, char sep, char lastSep);
/* Print out transMapGene.  Separate fields with sep. Follow last field with lastSep. */

#define transMapGeneTabOut(el,f) transMapGeneOutput(el,f,'\t','\n');
/* Print out transMapGene as a line in a tab-separated file. */

#define transMapGeneCommaOut(el,f) transMapGeneOutput(el,f,',',',');
/* Print out transMapGene as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

struct transMapGene *transMapGeneQuery(struct sqlConnection *geneConn,
                                       char *table, char *srcDb, char *srcId);
/* load a single transMapSrc object for an srcDb and srcId from a table,
 * or return NULL if not found */

#endif /* TRANSMAPGENE_H */

