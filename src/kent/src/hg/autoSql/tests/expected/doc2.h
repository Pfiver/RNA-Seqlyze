/* doc2.h was originally generated by the autoSql program, which also 
 * generated doc2.c and doc2.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef DOC2_H
#define DOC2_H

#define POINT_NUM_COLS 3

struct point
/* A three dimensional point */
    {
    float x;	/* Horizontal coordinate */
    float y;	/* Vertical coordinate */
    float z;	/* In/out of screen coordinate */
    };

struct point *pointLoad(char **row);
/* Load a point from row fetched with select * from point
 * from database.  Dispose of this with pointFree(). */

struct point *pointLoadAll(char *fileName);
/* Load all point from whitespace-separated file.
 * Dispose of this with pointFreeList(). */

struct point *pointLoadAllByChar(char *fileName, char chopper);
/* Load all point from chopper separated file.
 * Dispose of this with pointFreeList(). */

#define pointLoadAllByTab(a) pointLoadAllByChar(a, '\t');
/* Load all point from tab separated file.
 * Dispose of this with pointFreeList(). */

struct point *pointCommaIn(char **pS, struct point *ret);
/* Create a point out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new point */

void pointOutput(struct point *el, FILE *f, char sep, char lastSep);
/* Print out point.  Separate fields with sep. Follow last field with lastSep. */

#define pointTabOut(el,f) pointOutput(el,f,'\t','\n');
/* Print out point as a line in a tab-separated file. */

#define pointCommaOut(el,f) pointOutput(el,f,',',',');
/* Print out point as a comma separated list including final comma. */

void pointJsonOutput(struct point *el, FILE *f);
/* Print out point in JSON format. */

#define COLOR_NUM_COLS 3

struct color
/* A red/green/blue format color */
    {
    unsigned char red;	/* Red value 0-255 */
    unsigned char green;	/* Green value 0-255 */
    unsigned char blue;	/* Blue value 0-255 */
    };

struct color *colorLoad(char **row);
/* Load a color from row fetched with select * from color
 * from database.  Dispose of this with colorFree(). */

struct color *colorLoadAll(char *fileName);
/* Load all color from whitespace-separated file.
 * Dispose of this with colorFreeList(). */

struct color *colorLoadAllByChar(char *fileName, char chopper);
/* Load all color from chopper separated file.
 * Dispose of this with colorFreeList(). */

#define colorLoadAllByTab(a) colorLoadAllByChar(a, '\t');
/* Load all color from tab separated file.
 * Dispose of this with colorFreeList(). */

struct color *colorCommaIn(char **pS, struct color *ret);
/* Create a color out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new color */

void colorOutput(struct color *el, FILE *f, char sep, char lastSep);
/* Print out color.  Separate fields with sep. Follow last field with lastSep. */

#define colorTabOut(el,f) colorOutput(el,f,'\t','\n');
/* Print out color as a line in a tab-separated file. */

#define colorCommaOut(el,f) colorOutput(el,f,',',',');
/* Print out color as a comma separated list including final comma. */

void colorJsonOutput(struct color *el, FILE *f);
/* Print out color in JSON format. */

#define FACE_NUM_COLS 3

struct face
/* A face of a three dimensional solid */
    {
    struct face *next;  /* Next in singly linked list. */
    struct color color;	/* Color of this face */
    int pointCount;	/* Number of points in this polygon */
    unsigned *points;	/* Indices of points that make up face in polyhedron point array */
    };

struct face *faceLoad(char **row);
/* Load a face from row fetched with select * from face
 * from database.  Dispose of this with faceFree(). */

struct face *faceLoadAll(char *fileName);
/* Load all face from whitespace-separated file.
 * Dispose of this with faceFreeList(). */

struct face *faceLoadAllByChar(char *fileName, char chopper);
/* Load all face from chopper separated file.
 * Dispose of this with faceFreeList(). */

#define faceLoadAllByTab(a) faceLoadAllByChar(a, '\t');
/* Load all face from tab separated file.
 * Dispose of this with faceFreeList(). */

struct face *faceCommaIn(char **pS, struct face *ret);
/* Create a face out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new face */

void faceFree(struct face **pEl);
/* Free a single dynamically allocated face such as created
 * with faceLoad(). */

void faceFreeList(struct face **pList);
/* Free a list of dynamically allocated face's */

void faceOutput(struct face *el, FILE *f, char sep, char lastSep);
/* Print out face.  Separate fields with sep. Follow last field with lastSep. */

#define faceTabOut(el,f) faceOutput(el,f,'\t','\n');
/* Print out face as a line in a tab-separated file. */

#define faceCommaOut(el,f) faceOutput(el,f,',',',');
/* Print out face as a comma separated list including final comma. */

void faceJsonOutput(struct face *el, FILE *f);
/* Print out face in JSON format. */

#define POLYHEDRON_NUM_COLS 4

struct polyhedron
/* A solid three dimensional object */
    {
    struct polyhedron *next;  /* Next in singly linked list. */
    int faceCount;	/* Number of faces */
    struct face *faces;	/* List of faces */
    int pointCount;	/* Number of points */
    struct point *points;	/* Array of points */
    };

struct polyhedron *polyhedronLoad(char **row);
/* Load a polyhedron from row fetched with select * from polyhedron
 * from database.  Dispose of this with polyhedronFree(). */

struct polyhedron *polyhedronLoadAll(char *fileName);
/* Load all polyhedron from whitespace-separated file.
 * Dispose of this with polyhedronFreeList(). */

struct polyhedron *polyhedronLoadAllByChar(char *fileName, char chopper);
/* Load all polyhedron from chopper separated file.
 * Dispose of this with polyhedronFreeList(). */

#define polyhedronLoadAllByTab(a) polyhedronLoadAllByChar(a, '\t');
/* Load all polyhedron from tab separated file.
 * Dispose of this with polyhedronFreeList(). */

struct polyhedron *polyhedronCommaIn(char **pS, struct polyhedron *ret);
/* Create a polyhedron out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new polyhedron */

void polyhedronFree(struct polyhedron **pEl);
/* Free a single dynamically allocated polyhedron such as created
 * with polyhedronLoad(). */

void polyhedronFreeList(struct polyhedron **pList);
/* Free a list of dynamically allocated polyhedron's */

void polyhedronOutput(struct polyhedron *el, FILE *f, char sep, char lastSep);
/* Print out polyhedron.  Separate fields with sep. Follow last field with lastSep. */

#define polyhedronTabOut(el,f) polyhedronOutput(el,f,'\t','\n');
/* Print out polyhedron as a line in a tab-separated file. */

#define polyhedronCommaOut(el,f) polyhedronOutput(el,f,',',',');
/* Print out polyhedron as a comma separated list including final comma. */

void polyhedronJsonOutput(struct polyhedron *el, FILE *f);
/* Print out polyhedron in JSON format. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* DOC2_H */

