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
 * File:    algorithm.c
 *
 * This module provide the high-level algorithm for finding and
 * identifying all homologyteams.   The techniques and procedures used are
 * closely modeled after those in the paper
 *    "Identifying Conserved Gene Clusters in the Presence of
 *     Orthologous Groups", by Xin He and Michael H. Goldwasser,
 *     appearing in the Proceedings for RECOMB 2004.
 *
 ****************************************************************/
#include <stdlib.h>
#include "algorithm.h"
#include "filter.h"
#include "globaldata.h"
#include "list.h"
#include "types.h"

/***********************************************************
 * Global data used internally by the algorithm
 ***********************************************************/
int globalTime;
int* stamp;            // Array of timestamps
int* tempMark;         // Additional array for use within markCommonAlphabet
int* reconstructMark;  // Additional array used for reconstructing teams
Subchrom inProgressC; 
Subchrom inProgressD; 
Chromosome *identifyC; // Used to differentiate between the two originals

/***********************************************************
 * some utility functions
 ***********************************************************/

/*
 *  returns the distance between specified gene positions in subchrom.
 *  (assumes that pos1 <= pos2)
 */
/* [MP] : no need of this anymore
double dist(Subchrom* s, int pos1, int pos2) {
  return ( s->chrom->genes[pos2].position - s->chrom->genes[pos1].position );
}
*/


// start [MP] code here
/*
 *  returns the distance between specified gene positions in subchrom.
 *  (assumes that pos1 <= pos2)
 */
boole booldist(Subchrom* s, int pos1, int pos2) {
  switch (bparam) {
  case 0: return(s->chrom->genes[pos2].positionBeg - s->chrom->genes[pos1].positionEnd <= delta);
  case 1: return(s->chrom->genes[pos2].geneord - s->chrom->genes[pos1].geneord <= no_genes);
  case 2: return((s->chrom->genes[pos2].positionBeg - s->chrom->genes[pos1].positionEnd <= delta) &&
		 (s->chrom->genes[pos2].geneord - s->chrom->genes[pos1].geneord <= no_genes));
  default: return(0);
  }
}
// end [MP] code here

/*
 *  This checks whether the given subchromosomes are identical
 */
boole subEqual(Subchrom s1, Subchrom s2) {
  return(s1.chrom == s2.chrom && s1.startIndex==s2.startIndex && s1.endIndex==s2.endIndex);
}


/*
 *  This checks whether the given subchromosome is empty
 */
int subEmpty(Subchrom* sub) {
  return (sub->startIndex > sub->endIndex);
}


/*
 *  Assumes that excluded is a prefix of shrinking.
 *  This modifies 'shrinking' so as to exclude 'excluded'
 */
void subMinus(Subchrom* shrinking, Subchrom* excluded) {
  shrinking->startIndex = 1 + excluded->endIndex;
}


/*
 *  Allocates a new subchromosome which is set to represent the entire
 *  original chromosome.
 */
Subchrom subEntire(Chromosome* original) {
  Subchrom sub;
  sub.chrom = original;
  sub.startIndex = 0;
  sub.endIndex = original->numGenes - 1;
  return sub;
}


/*
 *  Allocates a new subchromosome which is set to represent the prefix
 *  of the original, up to and including the specified index.
 */
Subchrom subPrefix(Subchrom* s, int index) {
  Subchrom sub;
  sub = (*s);
  sub.endIndex = index;
  return sub;
}


/*
 * Used to explicitly free the global arrays created by direct calls
 * to initializeTime.   (the global arrays that were created for use
 * in a call to findTeams will be deallocated by that routine).
 */
void freeGlobal() {
  if (stamp) free(stamp);
  if (tempMark) free(tempMark);
  if (reconstructMark) free(reconstructMark);
}


/*
 * Initalizes the global arrays and resets the global time.
 */
void initializeTime()
{
  int c;
  freeGlobal();
  globalTime = 0;
  stamp    = (int*) malloc(num_families*sizeof(int));
  tempMark = (int*) malloc(num_families*sizeof(int));
  reconstructMark = (int*) malloc(num_families*sizeof(int));
  for (c=0; c<num_families; c++) {
    stamp[c] = tempMark[c] = reconstructMark[c] = 0;
  }
}

/*
 * Create a Delimiter structure, determining which source is C and which D.
 */
Delimiter delim(Subchrom *A, Subchrom *B)
{
  Delimiter del;

  if (A->chrom == identifyC) {
    del.subC = (*A);
    del.subD = (*B);
  } else {
    del.subC = (*B);
    del.subD = (*A);
  }

  inProgressC = del.subC;
  inProgressD = del.subD;

  return del;
}



/***********************************************************
 * Routines as outlined in paper
 ***********************************************************/

/*
 * Updates the global time, and then marks all famililies which are common
 * to the two subchromosomesw with the global time.
 *
 * The updated time is returned, for reference.
 */
int markCommonAlphabet(Subchrom A, Subchrom B) {
  int g;
  globalTime++;

  for (g = A.startIndex; g <= A.endIndex; g++)
    tempMark[A.chrom->genes[g].family] = globalTime;
  for (g = B.startIndex; g <= B.endIndex; g++)
    if (tempMark[B.chrom->genes[g].family] == globalTime)
      stamp[B.chrom->genes[g].family] = globalTime;
  return globalTime;
}


/*
 * Note, if two chromosomes have no common families, then the returned
 * subchromosome will be equivalent to S in entirety.
 */
Subchrom findFirstRun(Subchrom* S, int timestamp) {
  int endRun, nextGene;

  endRun = S->startIndex;
  while ((endRun < S->endIndex) && (stamp[S->chrom->genes[endRun].family] < timestamp))
    endRun++;
  nextGene = endRun+1;
  // before: while ((nextGene <= S->endIndex) && (dist(S,endRun,nextGene) <= delta)) {
  while ((nextGene <= S->endIndex) && booldist(S,endRun,nextGene)) { // [MP] version
    // check whether nextGene is in common alphabet,
    //   and thus whether to extend the run
    if (stamp[S->chrom->genes[nextGene].family] >= timestamp)
      endRun = nextGene;

    // in either case, continue advancing
    nextGene++;
  }

  return subPrefix(S,endRun);
}


/*
 * Recursive procedure to Find the homolgy teams of two subchromosomes.
 * Originally invoked via findTeams()
 */
void findTeamsRecurse(Subchrom A, Subchrom B) {

  int localTime = markCommonAlphabet(A,B);



  Subchrom A_first = findFirstRun(&A,localTime);  // the first delta-run of A
  subMinus(&A,&A_first);		// now A is the rest of the original A

  if (subEmpty(&A)) {

    filters[chosenFilter].reportTeam(delim(&A_first,&B));
    inProgressC.chrom = inProgressD.chrom = NULL;

  } else {            

    do {
      findTeamsRecurse(B,A_first);
      A_first = findFirstRun(&A,localTime);
      subMinus(&A,&A_first);
    } while (!subEmpty(&A));

    findTeamsRecurse(B,A_first);  // Note: fixes bug in RECOMB paper

  }
}


/*
 * The main algorithm that find homology teams
 */
void findTeams(Chromosome* origC, Chromosome* origD) {
  int localTime;
  Subchrom C_first, C_rest, D;
  identifyC = origC;

  // Initialize global data
  initializeTime();
  filters[chosenFilter].begin();

  //printf( "Start findTeams...\n" );
  // Decompose problem based on initial runs of C
  C_rest = subEntire(origC);   // convert orig chromosomes to
  D = subEntire(origD);        //   equivalent subchrom structures
  localTime = markCommonAlphabet(C_rest,D);
  //printf( "start\n" );

  do {

    C_first = findFirstRun(&C_rest,localTime);

    subMinus(&C_rest,&C_first);

    findTeamsRecurse(D,C_first);

  } while (!subEmpty(&C_rest));


  // deallocate global data
  filters[chosenFilter].end();
  freeGlobal();
}



/***********************************************************
 * Additional routine for reconstucting teams/witnesses
 ***********************************************************/

/*
 * This method is used as an aid to the reporting filters.
 *
 * Reconstructs the team common two the reported subchromsomes.  
 * (this assumes that the given subchrosomes have already been earlier   
 * identified as defining a true team)
 *
 * If selectedFamily >= 0, then will only return teams which contain that family.
 *
 * Returns: the number of families in the team (as opposed to genes)
 * Effect:  If contents!=NULL, will set it to reference a new array of
 *          ints, representing the families in the team.
 */
int reconstructTeam(Delimiter del, int selectedFamily, int** contents) {
  int i,family;
  int cardinality=0;
  Subchrom minLength;
  List *l;
  Node *node;

  if (!subEqual(del.subC,inProgressC) || !subEqual(del.subD,inProgressD)) {
    // presumed to be called AFTER completion of algorithm.
    // must recompute common alphabet
    markCommonAlphabet(del.subC,del.subD);
  } 

  if (contents)
    l = listNew();


  if (selectedFamily==-1 || stamp[selectedFamily] == globalTime) {
    minLength = ((( del.subC.endIndex - del.subC.startIndex) <
		  ( del.subD.endIndex - del.subD.startIndex)) ?
		 del.subC : del.subD);
    for (i=minLength.startIndex; i<=minLength.endIndex; i++) {
      family = minLength.chrom->genes[i].family;
      if (stamp[family] == globalTime) // its common
	if (reconstructMark[family] != globalTime) {
	  // its the first occurrence of this family
	  reconstructMark[family]=globalTime;
	  cardinality++;
	  if (contents) {
	    int* newI;
	    newI = (int*) malloc(sizeof(int));
	    *newI = family;
	    listAddBack(l,newI);
	  }
	}
    }
  }

  if (contents) {
    (*contents) = (int*) malloc(cardinality*sizeof(int));

    for (i=0,node=l->head->next; i<cardinality; i++) {
      int* tempI;
      tempI = (int*) node->item;
      (*contents)[i] = *tempI;
      free(tempI);
      node = node->next;
    }
    listFree(l);

  }

  return cardinality;
}


/*
 * Prints the witness.   This routine should only be called when it
 * is known that markCommonAlphabet was most recently called for this team.
 */
void printWitness(FILE* fp, Subchrom* A) {
  boole first=T;
  int i,family;
  for (i=A->startIndex; i<=A->endIndex; i++) {
    family = A->chrom->genes[i].family;
    if (stamp[family] == globalTime) {
      if (first)
	first = F;
      else
	fprintf(fp,",");
      // before:      fprintf(fp,"(%s:%0.2lf)",family_names[family],A->chrom->genes[i].position);
      fprintf(fp,"(%s:%0.2ld)",family_names[family],A->chrom->genes[i].positionBeg); // [MP] code
      //fprintf(fp,"(%s:%0.2ld:%0.2ld:%d)",family_names[family],A->chrom->genes[i].positionBeg,A->chrom->genes[i].positionEnd,A->chrom->genes[i].geneord); // [MP] code
    }
  }
  fprintf(fp,"\n");
}



