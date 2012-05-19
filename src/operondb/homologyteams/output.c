/****************************************************************
 *
 * Project: homologyteams
 * Authors: Michael H. Goldwasser (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
 * Version: 1.0 (May 2004)
 *
 * Copyright (C) 2004  Michael H. Goldwasser and Xin He
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *****************************************************************
 *
 * File:    output.c
 *
 * Support for outputting a desired team, in an appropriate format,
 * where a team is identified in delimited form based upon the two
 * identifying subchromosomes.  Currently, there exists support for
 * three formats:
 *     DELIM  - abbreviated output of original delimters
 *     FAMILIES   - output the homology families
 *     WITNESS - output the families, as well as the witnesses
 *
 ****************************************************************/
#include <string.h>
#include <stdlib.h>
#include "algorithm.h"
#include "output.h"
#include "types.h"


enum OutputType {FAMILIES=0, WITNESS, DELIM};
enum OutputType outputStyle = FAMILIES;   // Default
FILE* ofp;

typedef struct style_ {
  char* name;
  enum OutputType enumVal;
  char* description;
} style;

style outputStyles[] = {
  {"families",FAMILIES,"output for each team is list of homolgy families"},
  {"witness",WITNESS,"output for each team includes witnesses"},
  {"concise",DELIM,"output is concise summary of teams (to be later used with -T option)"}
};


/*
 *  Generates output for the team defined by the give delimiter.  
 *  The style of the output is based upon the current setting of
 *  outputStyle.
 *
 *  If not NULL, the second parameter is an array of the family IDs,
 *  with the third parameter the size of the array.
 */ 
void outputTeam(Delimiter delim, int** knownContents, int knownSize) {
  int cardinality,i;
  int* reconstruct;
  int** contents;


  switch (outputStyle) {
  case DELIM:
    fprintf(ofp,"%9d %9d %9d %9d\n",
	    delim.subC.startIndex,
	    delim.subC.endIndex,
	    delim.subD.startIndex,
	    delim.subD.endIndex);
    break;


  case WITNESS:
  case FAMILIES:
    if (knownContents) {
      contents = knownContents;
      cardinality = knownSize;
    } else {
      cardinality = reconstructTeam(delim,-1,&reconstruct);
      contents = &reconstruct;
    }

    if (cardinality>1) {

      for (i=0; i<cardinality-1; i++)
	fprintf(ofp,"%s:",family_names[(*contents)[i]]);
      fprintf(ofp,"%s\n",family_names[(*contents)[cardinality-1]]);

      if (outputStyle==WITNESS) {
	fprintf(ofp,"Witness in the 1st chromosome:\n   ");
	printWitness(ofp,&delim.subC);

	fprintf(ofp,"Witness in the 2nd chromosome:\n   ");
	printWitness(ofp,&delim.subD);

      }

      fprintf(ofp,"\n");
    }

    break;
  }

  fflush(ofp);
}



/*
 *  Sets the desired output style, based on the given 'key'.
 *
 *  returns failure if key not found
 */
boole chooseOutputStyle(char* key) {
  boole found = F;
  int entries, i;

  // loop through array of outputStyles, searching for key
  entries = sizeof(outputStyles)/sizeof(style);
  for (i=0; !found && i<entries; i++)
    if (strcmp(key,outputStyles[i].name)==0) {
      found = T;
      outputStyle = outputStyles[i].enumVal;
    }

  if (!found)
    fprintf(stderr,"Unrecognized output style (%s)\n",key);

  return found;
}



/*
 *  Prints the keywords for the existing output choices
 */
void printOutputStyles(char* preface) {
  int entries, i;

  entries = sizeof(outputStyles)/sizeof(style);

  printf("%s%-7s %s\n",preface,"style","description");
  printf("%s%-7s %s\n",preface,"-------","---------------------------");
  

  for (i=0; i<entries; i++) {
    printf("%s%-7s %s",preface,outputStyles[i].name,outputStyles[i].description);
    if (i==0)
      printf("  (DEFAULT)");
    printf("\n");
  }


}

