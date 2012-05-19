/* syntenySanger.h was originally generated by the autoSql program, which also 
 * generated syntenySanger.c and syntenySanger.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef SYNTENYSANGER_H
#define SYNTENYSANGER_H

struct syntenySanger
/* Sanger Mouse Synteny */
    {
    struct syntenySanger *next;  /* Next in singly linked list. */
    char *chrom;	/* Human Chrom */
    unsigned chromStart;	/* Start on Human */
    unsigned chromEnd;	/* End on Human */
    char *name;	/* Mouse Chromosome */
    unsigned score;	/* score always zero */
    char strand[2];	/* + direction matches - opposite */
    };

void syntenySangerStaticLoad(char **row, struct syntenySanger *ret);
/* Load a row from syntenySanger table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct syntenySanger *syntenySangerLoad(char **row);
/* Load a syntenySanger from row fetched with select * from syntenySanger
 * from database.  Dispose of this with syntenySangerFree(). */

struct syntenySanger *syntenySangerLoadAll(char *fileName);
/* Load all syntenySanger from a tab-separated file.
 * Dispose of this with syntenySangerFreeList(). */

struct syntenySanger *syntenySangerLoadWhere(struct sqlConnection *conn, char *table, char *where);
/* Load all syntenySanger from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with syntenySangerFreeList(). */

struct syntenySanger *syntenySangerCommaIn(char **pS, struct syntenySanger *ret);
/* Create a syntenySanger out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new syntenySanger */

void syntenySangerFree(struct syntenySanger **pEl);
/* Free a single dynamically allocated syntenySanger such as created
 * with syntenySangerLoad(). */

void syntenySangerFreeList(struct syntenySanger **pList);
/* Free a list of dynamically allocated syntenySanger's */

void syntenySangerOutput(struct syntenySanger *el, FILE *f, char sep, char lastSep);
/* Print out syntenySanger.  Separate fields with sep. Follow last field with lastSep. */

#define syntenySangerTabOut(el,f) syntenySangerOutput(el,f,'\t','\n');
/* Print out syntenySanger as a line in a tab-separated file. */

#define syntenySangerCommaOut(el,f) syntenySangerOutput(el,f,',',',');
/* Print out syntenySanger as a comma separated list including final comma. */

#endif /* SYNTENYSANGER_H */

