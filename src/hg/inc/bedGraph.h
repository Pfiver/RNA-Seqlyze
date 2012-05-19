/* bedGraph.h was originally generated by the autoSql program, which also 
 * generated bedGraph.c and bedGraph.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef BEDGRAPH_H
#define BEDGRAPH_H

#define BEDGRAPH_NUM_COLS 4

struct bedGraph
/* bed-like graphing data */
    {
    struct bedGraph *next;  /* Next in singly linked list. */
    char *chrom;	/* Chromosome or FPC contig */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    float dataValue;	/* data value for this range */
    };

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* BEDGRAPH_H */


