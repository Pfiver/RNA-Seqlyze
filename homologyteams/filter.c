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
 * File:    filter.c
 *
 * This file includes a variety of filters which control the precise
 * rules used in selecting which teams are included in the output.
 *
 * Additional such filters can be designed and included.  To do so,
 * the three appropriate functions must be written, and the filter
 * must be entered into the array of filters, with a unique identifier.
 *
 ****************************************************************/
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "algorithm.h"
#include "filter.h"
#include "types.h"
#include "output.h"



/****************************************************************
 *
 * The index of the chosen filter
 *
 ****************************************************************/
int chosenFilter;


/****************************************************************
 *
 * An optional parameter provided to the chosen filter
 *
 ****************************************************************/
char* filterParam;

/****************************************************************
 *
 * Some interpretations of the parameter
 *
 ****************************************************************/
int sizeThreshold;
int selectedFamily;


/****************************************************************
 *
 * Utilities
 *
 ****************************************************************/


/*
 * Find the number of a family in the family_names table
 * Return: -1 if not found; otherwise the number                   
 */
int findFamilyNumber( const char* name)
{
  int i;
  for (i=num_families-1; i>=0 && strcmp(family_names[i],name); i--);
  return i;
}



/****************************************************************
 *
 * This filter will report all teams.  It is the default.  (an attempt
 * is made to remove teams which are clearly trivial, but for more
 * control, should rely on one of the other filters.
 *
 ****************************************************************/
boole allCheck() {
  return T;
}

void allBegin() {

}

void allReportTeam(Delimiter delim) {
  if ((delim.subC.startIndex != delim.subC.endIndex) &&
      (delim.subD.startIndex != delim.subD.endIndex))
    outputTeam(delim,NULL,0);
}

void allEnd() {

}



/****************************************************************
 *
 * This filter will report only those teams with at least a given cardinality
 *
 ****************************************************************/
boole sizeCheck() {
  boole success = T;
  char* temp;

  if (filterParam) {
    sizeThreshold = (int) strtol(filterParam,&temp,10);
    success = (boole) !(sizeThreshold<1 || *temp!='\0');
    if (!success) {
      fprintf(stderr,"ERROR: parameter (%s) not recognized by filter %s.\n",
	      filterParam,filters[chosenFilter].name);
    }
  } else {
    fprintf(stderr,"ERROR: filter %s requires an integer parameter.\n",
	    filters[chosenFilter].name);
    success = F;
  }

  return success;
}

void sizeBegin() {

}

void sizeReportTeam(Delimiter delim) {
  int cardinality;
  int* contents;
  cardinality = reconstructTeam(delim,-1,&contents);

  if (cardinality>=sizeThreshold)
    outputTeam(delim,&contents,cardinality);

  free (contents);
}

void sizeEnd() {

}



/****************************************************************
 *
 * This filter will report only those teams which contain the
 * specified family.
 *
 ****************************************************************/
boole familyCheck() {
  // can check for the existence of the parameter, but not yet for
  // whether the given string is a true identifier because we have not
  // yet read the data from the file.  Will have to wait until
  // familyBegin to do this.


  boole success = T;
  if (!filterParam) {
    success = F;
    fprintf(stderr,"ERROR: filter %s requires that you specify a named family.\n",
	    filters[chosenFilter].name);
  }

  return (success);
}


void familyBegin() {

  if ((selectedFamily = findFamilyNumber(filterParam))==-1) {
    fprintf(stderr,"WARNING: specified  family name '%s' not recognized.\n",filterParam);
  }

}

void familyReportTeam(Delimiter delim) {
  int cardinality;
  int* contents;
  cardinality = reconstructTeam(delim,selectedFamily,&contents);

  if (cardinality>=-1)
    outputTeam(delim,&contents,cardinality);

  free (contents);
}

void familyEnd() {

}


/****************************************************************
 *
 * For each family, this filter will report the maximum cardinality
 * team which contains that particular family.
 *
 * Optionally, a single family can be specified, in which case results
 * will only be compiled for that single family.
 ****************************************************************/

typedef struct candidate_ {
  int cardinality;
  Delimiter delim;
} candidate;

candidate* best;

boole maxCheck() {
  return T;
}


void maxBegin() {
  int i;

  // check to see if optional parameter was specified
  selectedFamily = -1;  // default is to be interested in all families
  if (filterParam) {
    if ((selectedFamily = findFamilyNumber(filterParam))==-1) {
      fprintf(stderr,"WARNING: specified  family name '%s' not recognized. Will report results for each family.\n",filterParam);
    }
  }


  // set up table for compiling results
  best = (candidate*) malloc(num_families * sizeof(candidate));
  for (i=0; i<num_families; i++)
    best[i].cardinality = 0;
}


void maxReportTeam(Delimiter delim) {
  int cardinality,i;
  int* contents;
  cardinality = reconstructTeam(delim,selectedFamily,&contents);

  // record results if new best for each family in contents
  for (i=0; i<cardinality; i++)
    if (cardinality>best[contents[i]].cardinality){
      best[contents[i]].cardinality = cardinality;
      best[contents[i]].delim = delim;
    }

  free (contents);
}


void maxEnd() {
  if (selectedFamily != -1) {
    fprintf(ofp,"Of those teams containing family %s, the largest is\n",family_names[selectedFamily]);
    outputTeam(best[selectedFamily].delim,NULL,0);
  } else {
    int i;
    
    for (i=0; i<num_families; i++) {
      if (best[i].cardinality>1) {
	fprintf(ofp,"Of those teams containing %s, the largest has cardinality %d\n",family_names[i],best[i].cardinality);
	outputTeam(best[i].delim,NULL,0);
      }     
    }
  }

  free (best);
}





/****************************************************************
 *
 * The array of possible filters from which to choose
 *
 ****************************************************************/
filter filters[] = {
  {"all",  &allCheck,  &allBegin,  &allReportTeam,  &allEnd, "all teams",""},
  {"size", &sizeCheck, &sizeBegin, &sizeReportTeam, &sizeEnd, "teams with at least given cardinality","<val>"},
  {"family", &familyCheck, &familyBegin, &familyReportTeam, &familyEnd, "teams which contain the specified family","<family name>"},
  {"max", &maxCheck, &maxBegin, &maxReportTeam, &maxEnd, "for each family, the containing team with max cardinality",""},
  {"maxfamily", &maxCheck, &maxBegin, &maxReportTeam, &maxEnd, "for given family, the containing team with max cardinality","<family name>"}
};


/****************************************************************
 *
 * Some utility methods involving the array of filters
 *
 ****************************************************************/

/*
 *  Sets the chosenFilter variable, based on the given 'key'.  For some
 *  of the choices, param allows for specification of an optional
 *  parameter.
 *
 *  returns failure if key not found
 */
boole chooseFilter(char* key, char* param) {
  int entries, i;
  boole found = F;
  chosenFilter = 0;   // the default
  filterParam = param;

  // loop through array of filters, searching for key
  entries = sizeof(filters)/sizeof(filter);
  for (i=0; !found && i<entries; i++)
    if (strcmp(key,filters[i].name)==0) {
      found = T;
      chosenFilter = i;
    }

  if (found)
    found = filters[chosenFilter].checkParam();
  else
    fprintf(stderr,"Unrecognized filter (%s)\n",key);


  return found;
}



/*
 *  Prints the keywords and descriptions for the existing filters
 */
void printFilterChoices(char* preface) {

  int entries, i;

  entries = sizeof(filters)/sizeof(filter);

  printf("%s%-6s %-12s %s\n",preface,"rule","parameter","description");
  printf("%s%-6s %-12s %s\n",preface,"------","----------","---------------------------");
  

  for (i=0; i<entries; i++) {
    printf("%s%-6s %-12s %s",preface,filters[i].name,filters[i].paramdesc,filters[i].description);
    if (i==0)
      printf("  (DEFAULT)");
    printf("\n");
  }

}
