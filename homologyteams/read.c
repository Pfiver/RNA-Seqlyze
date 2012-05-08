/****************************************************************
 * This file has been modified from its original version -
 * homologyteams, version 1.0, provided by Michael H. Goldwasser 
 * (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
 * 
 * Copyright (C) 2008  Mihaela Pertea (mpertea@umiacs.umd.edu)
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
 * File:    read.c
 *
 * Support for reading files.
 *   -- file format for original definition of chromosomes pair
 *   -- file format for delimited team reports
 *
 ****************************************************************/
#include <stdlib.h>
#include <string.h>
#include "globaldata.h"
#include "types.h"
#include "list.h"
#include "read.h"

/*
 *  Used as comparison function for qsort of chromosomes.
 *
 *  Input: two genes
 *    Output: -1 if the location of 1st gene is less than 2nd gene
 *             0 if the location of 1st gene is equal to 2nd gene
 *             1 if the location of 1st gene is greater than 2nd gene
 */
int cmp_gene( const void* gene1, const void* gene2 )
{
  /* before :
  if ( ( (Gene*)gene1 )->position < ( (Gene*)gene2 )->position )
    return -1;
  if ( ( (Gene*)gene1 )->position == ( ( Gene*)gene2 )->position )
    return 0;
  else
    return 1;
  */
  // [MP] code start
  if ( ( (Gene*)gene1 )->positionBeg < ( (Gene*)gene2 )->positionBeg )
    return -1;
  if ( ( (Gene*)gene1 )->positionBeg == ( ( Gene*)gene2 )->positionBeg )
    return 0;
  else
    return 1;
  // [MP] code end
}



/******************************************************************************
 * This is a replacement for strsep which is not portable (missing on Solaris).
 */
static char* getToken(char** str, const char* delims)
{
  char* token;

  if (*str==NULL) {
    /* No more tokens */
    return NULL;
  }

  token=*str;
  while (**str!='\0') {
    if (strchr(delims,**str)!=NULL) {
      **str='\0';
      (*str)++;
      return token;
    }
    (*str)++;
  }
  /* There is no other token */
  *str=NULL;
  return token;
}




/*
 * read the original chromosome file, constructing both chromosomes, as well as the
 * family_names table.
 *
 * returns 0 if success, -1 if failure
 */
int readChromosomes( FILE* data_file, Chromosome* chrC, Chromosome* chrD, char*** familyArray )
{
  char buffer[MAX_LINE];
  char *line = buffer;
  char *delim = " \t\n";
  char *error;
  char *allPositions, *onePosition, *name, *keyword;
  char *temp; // [MP] code 
  List *listC, *listD, *listFamilyNames;
  long positionValBeg,positionValEnd; //[MP] code
  int geneOrder; //[MP] code
  int i;
  Node *node;

  // set these defaults in case of early exit due to error
  chrC->genes = NULL;
  chrD->genes = NULL;
  (*familyArray) = NULL;

  /* process the first line */
  if ( fgets(line, MAX_LINE, data_file ) == NULL )
    return -1;

  if (strlen(line)==MAX_LINE) {
    fprintf(stderr,"WARNING:  line %d of input file may exceed the maximum line length.",1);
  }

  /* the first word must be "FAMILY_ID", otherwise return -1 */
  keyword = getToken(&line, delim);
  if (keyword==NULL || strcmp( keyword, "FAMILY_ID" ) != 0 )
    return -1;


  name = getToken(&line,delim);
  chrC->name = strdup(name);

  name = getToken(&line,delim);
  chrD->name = strdup(name);



  /*
   * We will initially build chromosome as linked list of unordered genes.
   * Then will convert to an array and sort by position
   */
  listC = listNew();
  listD = listNew();
  listFamilyNames = listNew();



  num_families = 0;           // start from the first data line (2nd line of the data_file)
  while ( !feof( data_file ) ) {

    /* read one line at a time*/
    line = buffer;
    if ( fgets( line, MAX_LINE, data_file ) == NULL )
      continue;

    if (strlen(line)==MAX_LINE) {
      fprintf(stderr,"WARNING:  line %d of input file may exceed the maximum line length.",1+num_families);
    }

    /* tokenize the line */
    name = getToken(&line, delim);
    if ( *name == '\0' || !strcmp( name, "\n" ) )
      continue;


    /* record homology family name */
    listAddBack(listFamilyNames,strdup(name));

    //fprintf(stderr,"name=%s\n",name);

    /* process genes in C */
    allPositions = getToken( &line, delim );
    if ( strcmp( allPositions, "null" ) ) {

      while ((onePosition=getToken(&allPositions, ":")) != NULL) {
	temp=getToken(&onePosition,"-"); // start [MP] code

	if(temp !=NULL) {
	  positionValBeg =  strtol(temp, &error,10);

	  if (error == temp) {
	    // ERROR occurred
	    fprintf(stderr,"Error at line %d of input file.",(num_families+1));
	    fprintf(stderr,"Specified position (%s) is invalid.\n",temp);
	    fprintf(stderr,"Beg onePosition=%s\n",onePosition);
	  } else {
	    temp=getToken(&onePosition,"-"); 
	    if(temp !=NULL) {
	      positionValEnd =  strtol(temp, &error,10);

	      if (error == temp) {
		// ERROR occurred
		fprintf(stderr,"Error at line %d of input file.",(num_families+1));
		fprintf(stderr,"Specified position (%s) is invalid.\n",temp);
		fprintf(stderr,"End onePosition=%s\n",onePosition);
	      } else {
		temp=getToken(&onePosition,"-");
		if(temp != NULL) {
		  geneOrder=atoi(temp);

		  if(geneOrder<0) {
		    fprintf(stderr,"Error at line %d of input file.",(num_families+1));
		    fprintf(stderr,"Specified order (%s) is invalid.\n",temp);
		    fprintf(stderr,"Order onePosition=%s\n",onePosition);
		  }
		  else {
		    Gene *newGene = (Gene*) malloc(sizeof(Gene));
		    newGene->family = num_families;
		    newGene->positionBeg = positionValBeg;
		    newGene->positionEnd = positionValEnd;
		    newGene->geneord = geneOrder;

		    //fprintf(stderr,"1 beg=%ld end=%ld order=%d\n",positionValBeg,positionValEnd,geneOrder);

		    listAddBack(listC,newGene);
		  }
		}
	      }
	    }
	  } // end [MP] code here
        }
      }
    }


    /* process genes in D */
    allPositions = getToken( &line, delim );
    if ( strcmp( allPositions, "null" ) ) {

      while ((onePosition=getToken(&allPositions, ":")) != NULL) {
	temp=getToken(&onePosition,"-"); // start [MP] code
	if(temp != NULL) {
	  positionValBeg =  strtol(temp, &error,10);
	  if (error == temp) {
	    // ERROR occurred
	    fprintf(stderr,"Error at line %d of input file.",(num_families+1));
	    fprintf(stderr,"Specified position (%s) is invalid.\n",temp);
	    fprintf(stderr,"Beg onePosition=%s\n",onePosition);
	  } else {
	    temp=getToken(&onePosition,"-"); // start [MP] code
	    if(temp !=NULL) {
	      positionValEnd =  strtol(temp, &error,10);
	      if (error == temp) {
		// ERROR occurred
		fprintf(stderr,"Error at line %d of input file.",(num_families+1));
		fprintf(stderr,"Specified position (%s) is invalid.\n",temp);
		fprintf(stderr,"End onePosition=%s\n",onePosition);
	      } else {
		temp=getToken(&onePosition,"-");
		if(temp != NULL) {
		  geneOrder=atoi(temp);
		  if(geneOrder<0) {
		    fprintf(stderr,"Error at line %d of input file.",(num_families+1));
		    fprintf(stderr,"Specified order (%s) is invalid.\n",temp);
		    fprintf(stderr,"Order onePosition=%s\n",onePosition);
		  } else {
		    Gene *newGene = (Gene*) malloc(sizeof(Gene));
		    newGene->family = num_families;
		    newGene->positionBeg = positionValBeg;
		    newGene->positionEnd = positionValEnd;
		    newGene->geneord = geneOrder;

		    //fprintf(stderr,"2 beg=%ld end=%ld order=%d\n",positionValBeg,positionValEnd,geneOrder);

		    listAddBack(listD,newGene);
		  }
		}
	      } // end [MP] code
	    }
	  }
	}
      }
    }

    /* finished processing the current line */
    num_families++;
  }


  /* convert the various lists to arrays */

  (*familyArray) = (char**) malloc(sizeof(char*)*num_families);
  for (i=0,node=listFamilyNames->head->next; i<num_families; i++) {
    (*familyArray)[i] = (char*) node->item;
    node = node->next;
  }
  listFree(listFamilyNames);

  chrC->numGenes = listSize(listC);
  chrC->genes = (Gene*) malloc(chrC->numGenes*sizeof(Gene));
  for (i=0,node=listC->head->next; i<chrC->numGenes; i++) {
    Gene *temp = (Gene*) node->item;
    chrC->genes[i] = *temp;
    free(temp);
    node = node->next;
  }
  listFree(listC);

  qsort(chrC->genes, chrC->numGenes, sizeof(Gene), cmp_gene );

  chrD->numGenes = listSize(listD);
  chrD->genes = (Gene*) malloc(chrD->numGenes*sizeof(Gene));
  for (i=0,node=listD->head->next; i<chrD->numGenes; i++) {
    Gene *temp = (Gene*) node->item;
    chrD->genes[i] = *temp;
    free(temp);
    node = node->next;
  }
  listFree(listD);

  qsort(chrD->genes, chrD->numGenes, sizeof(Gene), cmp_gene );

  return 0;
}



/*
 * Reads next delimiter from the given file, and sets the given delim
 * appropriately.
 *
 * Returns true if success, false if failure
 */
boole readNextDelim(FILE* delimFile, Delimiter *delim) {
  boole success =
    (fscanf(delimFile,"%d %d %d %d\n",
            &delim->subC.startIndex,&delim->subC.endIndex,
            &delim->subD.startIndex,&delim->subD.endIndex) == 4);

  if (!success && !feof(delimFile))
    fprintf(stderr,"ERROR reading concise team file.\n");

  return success;
}
