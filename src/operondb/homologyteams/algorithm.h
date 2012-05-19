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
 * File:    algorithm.h
 *
 * This module provides the high-level findTeams algorithm,  as
 * described in the   
 *    "Identifying Conserved Gene Clusters in the Presence of
 *     Orthologous Groups", by Xin He and Michael H. Goldwasser,
 *     appearing in the Proceedings for RECOMB 2004.
 *
 * In addition, several related utlities are provided for use in the
 * recording module in reconstructing the actual teams which are
 * identified according to the defining subchromosomes.
 *
 ****************************************************************/
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <stdio.h>
#include "types.h"

/*
 * The main algorithm that find Family teams
 */
void findTeams(Chromosome* C, Chromosome* D);


/*
 * Initalizes the global arrays and resets the global time.  This does
 * not need to be explicitly called if using the findTeams routine.
 * However, it should be invoked in advance if the markCommonAlphabet
 * routine will be used independently.
 */
void initializeTime();


/*
 * This method is used as an aid to the reporting filters.
 *
 * Reconstructs the team common two the reported subchromsomes.  
 * (this assumes that the given subchrosomes have already been earlier   
 * identified as defining a true team)
 *
 * If selectedFamily >= 0, then will only return teams which contain that family.
 *
 * Returns: the number of Familys in the team (as opposed to genes)
 * Effect:  If contents!=NULL, will set it to reference a new array of
 *          ints, representing the Familys in the team.
 */
int reconstructTeam(Delimiter del, int selectedFamily, int** contents);

/*
 * Prints the witness.   This routine should only be called when it
 * is known that markCommonAlphabet was most recently called for this team.
 */
void printWitness(FILE* fp, Subchrom* A);

/*
 * Used to explicitly free the global arrays created by direct calls
 * to initializeTime.   (the global arrays that were created for use
 * in a call to findTeams will be deallocated by that routine).
 */
void freeGlobal();


#endif
